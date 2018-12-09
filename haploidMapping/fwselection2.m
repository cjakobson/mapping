%%% Function for using forward selection to make a genotype to phenotype
%%% map, which is then be fine mapped via QTN score algorithm
%%% Richard She
%%% December 07 2016
function [] = fwselection2(masterPath,savePath,pcutoff,pcutoff2,p)

%%% All input variables from server are characters as default - but we want
%%% p and pcutoff to be type (double).
if ischar(p)
    p = str2num(p);
end
if ischar(pcutoff)
    pcutoff = str2num(pcutoff);
end
%%% Make sure all pathnames end with a slash
if masterPath(length(masterPath)) ~= '/'
    masterPath = [masterPath,'/'];
end
if savePath(length(savePath)) ~= '/'
    savePath = [savePath,'/'];
end

%%% Load in depedencies.
cd(masterPath);
load([masterPath,'nutrientFilename.mat']);
load([masterPath,'nutrientTrait.mat']);
load([masterPath,'recursiveBootstrapCutoff2.mat']);
%load([masterPath,'phasedGenotype_preEmptyWells.mat']);
load([masterPath,'phasedGenotype.mat']);
load([masterPath,'variantPos.mat']);

%%% Check if filename and trait are the same length
if length(filename) ~= length(trait)
    fprintf('Filename.mat and trait.mat are not the same length \n');
end

% Create "edge well" pseudo-genototype.
col1 = ones(12,1);
col2 = [1;0;0;0;0;0;0;0;0;0;0;1];
plate = [col1;repmat(col2,6,1);col1];
edges = repmat(plate,12,1);

%%% Create pseudo genotypes for each of 12 individual plates
plateNorm = zeros(96*12,12);
for i = 1:12
    plateNorm(96*(i-1) + 1:96*i,i) = 1;
end

%%% Add these pseudogenotypes to the real genotype so we can regress
%%% against them
phasedGenotype2 = [edges,plateNorm,phasedGenotype];

% Exclude wells that don't meet coverage threshold
variantCalls = sum(abs(phasedGenotype2),1);
s = sum(abs(phasedGenotype2),2);

%%% Cutoff for number of variant calls needed to include a segregant.
rowCutoff = 5500;
columnCutoff = 0;

traitRow = trait{p};

x=phasedGenotype2(s>rowCutoff & traitRow ~= min(traitRow),variantCalls>columnCutoff);
% Remove edge count and plateNorm from variant calls
variantCalls(1:13)=[];
y=traitRow(s>rowCutoff & traitRow ~= min(traitRow));

%remove missing data
x=x(y~=0,:);
y=y(y~=0);


tic

% regress with forward selection
[b_fwselection,se,pval,inmodel,stats,nextstep,history] = stepwisefit(x,y,'penter',pcutoff,'display','on','maxiter',round(length(y)/6));
dev_fwselection = 1-stats.SSresid/stats.SStotal;
dof_fwselection = stats.df0;
bPos = find(inmodel);

toc

stats = rmfield(stats,'xr');
stats = rmfield(stats,'covb');

stepHistory = zeros(length(bPos),1);
for i = 1:length(bPos)
    stepHistory(i) = min(find(history.in(:,bPos(i))==1));
end

% Create one unified data structure to save all relevant workspace
% variables
b.b_fwselection = b_fwselection;
b.bPos = bPos;
b.dev_fwselection = dev_fwselection;
b.stats = stats;
b.inmodel = inmodel;
b.y = y;
b.x = x;
b.history = stepHistory;
prefix = [savePath,filename{p}];
save([prefix,'_cutoff_',num2str(pcutoff),'_b_2.mat'],'b','-v7.3'); 

% Weighted least squares regression;
% stdReps = error{p};
% stdReps = stdReps(s>rowCutoff & traitRow' ~= -1);
% regressionWeight = zeros(length(stdReps),1);
% for i = 1:length(stdReps)
%     regressionWeight(i) = 1/max(median(stdReps),stdReps(i));
% end
% mdl = fitlm(x(:,bPos),y,'Weights',regressionWeight);
% tempMdlResi = table2array(mdl.Residuals);
% mdlResiduals = tempMdlResi(:,1);

% Find 1D lod score for regressing each predictor one at a time
lod1D = zeros(length(x),1);
for i = 1:length(x)
    r = corr(x(:,i),y);
    lod1D(i) = -length(y)*log(1-r^2)/(2*log(10));
end

% Split variantPos into chromosome and position variables so that we can
% scan 10kb on either side of each causal variant
allPos = zeros(length(variantPos),1);
allChrom = cell(length(variantPos),1);
for i = 1:length(variantPos)
    tempVar = strsplit(variantPos{i},':');
    allChrom{i} = tempVar{1};
    allPos(i) = str2num(tempVar{2});
end

nVariants = sum(stats.PVAL>pcutoff2);
bPos_filtered = find(-log10(stats.PVAL) > pcutoff2);
[~,~,r_filtered] = regress(y,[ones(length(y),1),x(:,bPos_filtered)]);

tic

startIndex = sum(bPos_filtered <= 13) + 1; % First 13 bPos rows are edge normalization and plate normalization.
for posIndex = startIndex:length(bPos_filtered)
    posIndex
    length(bPos_filtered)
    pos1 = bPos_filtered(posIndex);
    tempVar = strsplit(variantPos{pos1-13},':');
    chrom = tempVar{1};
    position = str2num(tempVar{2});
    % Scan all positions up to 10 kb in either direction for possible
    % causal variants that fit better by the allele swap criteria.
    lower = min(find(strcmp(chrom,allChrom) & allPos > position - 10000)) + 13;
    upper = max(find(strcmp(chrom,allChrom) & allPos < position + 10000)) + 13;
    
    
    for i = lower:upper
        for j = lower:upper
            
            [p_anova{posIndex}(i-lower+1,j-lower+1),pairwise_p{posIndex}{i-lower+1}{j-lower+1},...
                ph1{posIndex}(i-lower+1,j-lower+1),ph2{posIndex}(i-lower+1,j-lower+1)] = ...
                fineMappingLod_multiSite_anova(i,j,lower,upper,pos1,x,y,b_fwselection,r_filtered);
            
        end
    end
    
end

toc

b.lod1D = lod1D;
b.p_anova = p_anova;
b.pairwise_p = pairwise_p;
b.ph1 = ph1;
b.ph2 = ph2;
b.cutoff2 = cutoff2(p);

save([prefix,'_cutoff_',num2str(pcutoff),'_b_2.mat'],'b','-v7.3');


