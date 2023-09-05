%simulate phenotype then perform regression
%do N replicates and output performance

function [] = linearMixedModelRadRemapPerm(traitIdx,nPerms)

%Pathname='/Users/cjakobson/Dropbox/JaroszLab/yeastCrossNutrientScreen/crosscrossdata/regression/modelSelection/incremental'

% if Pathname(length(Pathname)) ~= '/'
%     Pathname = [Pathname,'/'];
% end

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

for i=1:nPerms
tic
rng(i)
scramble=phenotypes(randperm(length(phenotypes)));
[b_fwselection,se,pval,inmodel,stats,nextstep,history] = stepwisefit(modelGenotypes,scramble,'penter',10^-3,'display','off');
dev_fwselection = 1-stats.SSresid/stats.SStotal;
dof_fwselection = stats.df0;
bPos = find(inmodel);
dof = length(bPos);
pValues{i} = -log10(stats.PVAL(bPos));
[pValues{i},sortIndex] = sort(pValues{i},'descend');
bPos = bPos(sortIndex);
toc     %this fit takes about 25min on sherlock
end



% Remove variables that aren't needed that would clog up HD space for when
% we save
clear genotypes;
clear stats; clear se; clear pval; clear domB;
clear inmodel; clear inmodel2; clear inmodel3; clear domSubset; clear newResidual; 
clear history; clear phasedVLCgenotype; clear modelGenotypes;
clear secondOrderGenotype;

% Save all the variables
save(['radRemap/' filename{traitIdx} '_perm.mat']);

end



