%simulate phenotype then perform regression
%do N replicates and output performance

function [] = linearMixedModelRadRemap(traitIdx,doHets,doFineMapping)


if ischar(traitIdx)
    traitIdx=str2num(traitIdx);
end

traitIdx

%cd(Pathname);
mkdir('radRemap');


load('phasedVLCgenotype.mat')
genotypes=phasedVLCgenotype;
clear phasedVLCgenotype

load('radRemapTrait.mat')
load('radRemapFilename.mat')

phenotypes=trait{traitIdx};
phenotypes(isnan(phenotypes))=0;

nPlates=length(phenotypes)/384;
nStrains=length(phenotypes);

[~, nCols]=size(genotypes);
nLoci=nCols/4;

%truncate to measured strains
genotypes=genotypes(1:nStrains,:);

%zero out missing growth measurements and missing genotypes
vNoSpot=phenotypes==min(phenotypes);
vNoGenotype=sum(genotypes==0,2)==nCols;

phenotypes(vNoSpot)=0;
phenotypes(vNoGenotype)=0;
genotypes(vNoSpot,:)=zeros(sum(vNoSpot),nCols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MODEL A
%ignore hets; homoRM gets 1; homoYJM gets -1
%final: only nStrains columns
[~,temp]=size(genotypes);
nLoci=temp/4;    

modelGenotypes=zeros(nStrains,nLoci);
for i=1:nStrains
    vGenotype=genotypes(i,1:nLoci)-genotypes(i,(nLoci+1):(2*nLoci));
    modelGenotypes(i,:)=vGenotype;
end

clear genotypes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pseudogenotypes=zeros(nStrains,nPlates+1);

%plates
for i=1:nPlates
    
    vPlate=zeros(nStrains,1);
    vPlate(((i-1)*384+1):(384*i),1)=ones(384,1);
    pseudogenotypes(:,i)=vPlate;
    
end

%edges
%top
vEdge=zeros(384,1);
vEdge(1:24)=1;
%bottom
vEdge(361:384)=1;
%sides
for i=2:15
    vEdge(24*(i-1)+1)=1;
    vEdge(24*i)=1;
end

for i=1:nPlates
    pseudogenotypes(((i-1)*384+1):(384*i),nPlates+1)=vEdge;
end


modelGenotypes=[pseudogenotypes modelGenotypes];

size(modelGenotypes)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
[b_fwselection,se,pval,inmodel,stats,nextstep,history] = stepwisefit(modelGenotypes,phenotypes,'penter',10^-3,'display','on');
dev_fwselection = 1-stats.SSresid/stats.SStotal;
dof_fwselection = stats.df0;
bPos = find(inmodel);
dof = length(bPos);
pValues = -log10(stats.PVAL(bPos));
[pValues,sortIndex] = sort(pValues,'descend');
bPos = bPos(sortIndex);
toc     %this fit takes about 25min on sherlock

%now add in het terms and refit

if doHets

load('phasedVLCgenotype.mat')
genotypes=phasedVLCgenotype;
clear phasedVLCgenotype

nPlates=length(phenotypes)/384;
nStrains=length(phenotypes);

[~, nCols]=size(genotypes);
nLoci=nCols/4;

%truncate to measured strains
genotypes=genotypes(1:nStrains,:);

%zero out missing growth measurements and missing genotypes
vNoSpot=phenotypes==min(phenotypes);
vNoGenotype=sum(genotypes==0,2)==nCols;

genotypes(vNoSpot,:)=zeros(sum(vNoSpot),nCols);


nStrains=length(phenotypes);

%MODEL C
%homoRM gets 1; homoYJM gets -1
%hets all get 1
%final: 2*nLoci columns
[~,temp]=size(genotypes);
nLoci=temp/4;    

modelGenotypes=zeros(nStrains,2*nLoci);
for i=1:nStrains
    vGenotype=[genotypes(i,1:nLoci)-genotypes(i,(nLoci+1):(2*nLoci)) genotypes(i,(2*nLoci+1):(3*nLoci))+genotypes(i,(3*nLoci+1):(4*nLoci))];
    modelGenotypes(i,:)=vGenotype;
end

clear genotypes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pseudogenotypes=zeros(nStrains,nPlates+1);

%plates
for i=1:nPlates
    
    vPlate=zeros(nStrains,1);
    vPlate(((i-1)*384+1):(384*i),1)=ones(384,1);
    pseudogenotypes(:,i)=vPlate;
    
end

%edges
%top
vEdge=zeros(384,1);
vEdge(1:24)=1;
%bottom
vEdge(361:384)=1;
%sides
for i=2:15
    vEdge(24*(i-1)+1)=1;
    vEdge(24*i)=1;
end

for i=1:nPlates
    pseudogenotypes(((i-1)*384+1):(384*i),nPlates+1)=vEdge;
end

modelGenotypes=[pseudogenotypes modelGenotypes];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tic
inmodel=[inmodel,zeros(1,nLoci)];   %second genotype mat has nLoci more columns
[b_fwselection,se,pval,inmodel,stats,nextstep,history] = stepwisefit(modelGenotypes,phenotypes,'penter',10^-3,'inmodel',logical(inmodel),'display','off');
dev_fwselection = 1-stats.SSresid/stats.SStotal;
dof_fwselection = stats.df0;
bPos = find(inmodel);
dof = length(bPos);
pValues = -log10(stats.PVAL(bPos));
[pValues,sortIndex] = sort(pValues,'descend');
bPos = bPos(sortIndex);
toc     %this fit takes about 5 min
end
%now do fine mapping

tic

if doFineMapping
    %neglect geometric factors (plates, edges) now and merge back later
    posToMap=bPos(find((bPos>(nPlates+1)).*(pValues>3)'));
    
    %remove those too close to the end (can't map)
    posToMap=posToMap(posToMap<(12054*2-10+(nPlates+1)));
    posToMap=posToMap(posToMap>(10+(nPlates+1)));

    %calculate residuals for fine mapping
    [~,~,r] = regress(phenotypes,[ones(length(phenotypes),1),modelGenotypes(:,bPos)]);
    
    ph2=cell(length(posToMap),1);
    
    for k=1:length(posToMap)

        position1=posToMap(k);
        upper=position1+10;
        lower=position1-10;

        %this code from RS routine
        for i = lower:upper
            for j = lower:upper

                %in RS code, x is genotypes and y is phenotypes
                [ph2{k}(i-lower+1,j-lower+1)] = ...
                    fineMappingLod_multiSite_anova(i,j,position1,modelGenotypes,...
                    b_fwselection,r);

            end
        end

    end

toc
%interpret fine mapping 

candidates=cell(length(posToMap),1);

for i=1:length(ph2)
    
    [~,candidates{i}]=qtnScore(ph2{i});
    
end

vResolved=0;
for i=1:length(candidates)
    vResolved(i)=length(candidates{i})==1;
end

fracResolved=sum(vResolved)/length(vResolved);
end


%%% Calculate percentage of variance explained by each predictor in
%%% the model
sumR = zeros(length(bPos),1);
varianceExplained = zeros(length(bPos),1);
for i = 1:length(bPos)
    newResidual = stats.yr + b_fwselection(bPos(i))*modelGenotypes(:,bPos(i));
    sumR(i) = sum(newResidual.^2) - stats.SSresid;
end
for i = 1:length(bPos)
    varianceExplained(i) = sumR(i)/sum(sumR)*dev_fwselection;
end




% Remove variables that aren't needed that would clog up HD space for when
% we save
clear genotypes;
clear stats; clear se; clear pval; clear domB;
clear inmodel; clear inmodel2; clear inmodel3; clear domSubset; clear newResidual; 
clear history; clear phasedVLCgenotype; clear modelGenotypes;
clear secondOrderGenotype;

% Save all the variables
save(['radRemap/' filename{traitIdx} '.mat']);

end



