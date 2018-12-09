%%% Script for parsing outputs from function "fwselection.m" (b variable)
% Identifies and annotates causal variants

clear all; close all; clc;

%%%%% CHANGE THIS TO THE MASTER PATH WHERE .mat FILES ARE KEPT %%%%%%%%%%%%
pathname = '/Users/cjakobson/Dropbox/JaroszLab/yeastCrossNutrientScreen/popGenRevisions/experiments/phenotyping/mapping/';
if pathname(length(pathname)) ~= '/'
    pathname = [pathname,'/'];
end
load([pathname,'nutrientTrait.mat']) %%%% Your Raw 1x1152 trait values
load([pathname,'variantPos.mat']);
load([pathname,'recursiveBootstrapCutoff2.mat'])
%%% Load vcf file that matches all variants in variantPos. 
[vcf,header,af,filler] = read_vcf([pathname,'variantPos080316_trimmedGenotype.eff.vcf']);

% Load SGD feature in features file
annotationFolder = [pathname,'annotateSGD/'];
featuresFile = fopen([annotationFolder,'SGD_features.tab']);
features = textscan(featuresFile,repmat('%s',1,16),'delimiter',char(9));
fclose(featuresFile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(pathname);
files = dir('*_b_2.mat');
filename = cell(length(files),1);

%%%% Parse filenames for useful identifiers. I typically include the time
%%%% and drug condition in the filename, but if these are not included,
%%%% we'll make a dummy variable.
for i = 1:length(files)
    filename{i} = files(i).name;
    temp=filename{i};
    if ~isempty(regexp(temp,'-'));
        temp = strsplit(filename{i},'-');
        time{i} = temp{1};
        if length(temp) > 1 & ~isempty(regexp(temp{2},'_'))
            temp = strsplit(temp{2},'_');
            condition{i} = temp{1};
        else
            condition{i} = ['dummyCondition',num2str(i)];
        end
    else
        time{i} = ['dummyTime',num2str(i)];
        condition{i} = ['dummyCondition',num2str(i)];
    end
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

%%% Extract SnpEff annotation line
annotation = cell(length(vcf{1}),1);
subannotation = cell(length(vcf{1}),1);
classifier = cell(length(vcf{1}),1);
mutation = cell(length(vcf{1}),1);
ORFannotation = cell(length(vcf{1}),1);
ORF = cell(length(vcf{1}),1);

for i = 1:length(vcf{1})
    temp = strsplit(vcf{8}{i},'ANN=');
    annotation{i} = temp{2}(4:length(temp{2}));
    if isempty(strfind(annotation{i},','))
        subannotation{i} = annotation{i};
        eff_fields = strsplit(subannotation{i},'|');
        classifier{i} = eff_fields{2};
        mutation{i} = [];
        for k = 1:length(eff_fields)
            if ~isempty(strfind(eff_fields{k},'c.'))
                mutation{i} = [mutation{i},eff_fields{k},'|'];
            end
            if ~isempty(strfind(eff_fields{k},'p.'))
                mutation{i} = [mutation{i},eff_fields{k},'|'];
            end
        end
        
        genestarts = regexp(subannotation{i},'Y[A-Z][LR][0-9][0-9][0-9][CW]');
        geneEnds = strfind(subannotation{i},'|');
        if ~isempty(genestarts)
            geneEnd = geneEnds(sum(geneEnds<genestarts(1))+1);
            geneList = subannotation{i}(genestarts(1):geneEnd-1);
        else
            tempStart = 3;
            geneList = subannotation{i}(geneEnds(tempStart)+1:geneEnds(tempStart+1)-1);
            while isempty(geneList) && tempStart < length(geneEnds) - 2
                tempStart = tempStart + 1;
                geneList = subannotation{i}(geneEnds(tempStart)+1:geneEnds(tempStart)-1);
            end
        end
        
        ORFannotation{i} = [];
        query = geneList;
        genestarts = regexp(query,'Y[A-Z][LR][0-9][0-9][0-9][CW]');
        if length(genestarts) > 1
            query = query(1:genestarts(2)-2);
        end
        index = find(strcmp(query,features{4}));
        if isempty(index)
            index = find(strcmp(query,features{5}));
        end
        if ~isempty(index)
            ORFannotation{i} = [ORFannotation{i},features{4}{index},':',...
                features{5}{index},':',features{16}{index}];
        else
            ORFannotation{i} = query;
        end
    else
        subannotation{i} = strsplit(annotation{i},',');
        classifier{i} = [];
        mutation{i} = [];
        for j = 1:length(subannotation{i})
            eff_fields = strsplit(subannotation{i}{j},'|');
            classifier{i} = [classifier{i},eff_fields{2},'|'];
            
            %%% KEEP ORF from the first subannotation
            if j == 1
                genestarts = regexp(subannotation{i}{j},'Y[A-Z][LR][0-9][0-9][0-9][CW]');
                geneEnds = strfind(subannotation{i}{j},'|');
                if ~isempty(genestarts)
                    geneEnd = geneEnds(sum(geneEnds<genestarts(1))+1);
                    geneList = subannotation{i}{j}(genestarts(1):geneEnd-1);
                else
                    tempStart = 3;
                    geneList = subannotation{i}{j}(geneEnds(tempStart)+1:geneEnds(tempStart+1)-1);
                    while isempty(geneList) && tempStart < length(geneEnds) - 2
                        tempStart = tempStart + 1;
                        geneList = subannotation{i}{j}(geneEnds(tempStart)+1:geneEnds(tempStart)-1);
                    end
                end
            end
            
            for k = 1:length(eff_fields)
                if ~isempty(strfind(eff_fields{k},'c.'))
                    mutation{i} = [mutation{i},eff_fields{k},'|'];
                end
                if ~isempty(strfind(eff_fields{k},'p.'))
                    mutation{i} = [mutation{i},eff_fields{k},'|'];
                end
            end
        end
        ORFannotation{i} = [];
        query = geneList;
        genestarts = regexp(query,'Y[A-Z][LR][0-9][0-9][0-9][CW]');
        if length(genestarts) > 1
            query = query(1:genestarts(2)-2);
        end
        index = find(strcmp(query,features{4}));
        if isempty(index)
            index = find(strcmp(query,features{5}));
        end
        if ~isempty(index)
            ORFannotation{i} = [ORFannotation{i},features{4}{index},':',...
                features{5}{index},':',features{16}{index}];
        else
            ORFannotation{i} = query;
        end
    end
end
for i = 1:length(ORF)
    if ~isempty(ORFannotation{i});
        temp = strsplit(ORFannotation{i},':');
        if iscell(temp)
            ORF{i} = temp{1};
        else
            ORF{i} = temp;
        end
    end
end

%%
%%% Collect metrics from saved b variable (loop through all 53 files) and
%%% write relevant metrics to a causal variants file

% Load in workspace from previous block.
pathname = '/Users/cjakobson/Dropbox/JaroszLab/yeastCrossNutrientScreen/popGenRevisions/experiments/phenotyping/mapping/';
cd(pathname);
files = dir('*0.001_b_2.mat');
filename = cell(length(files),1);
for i = 1:length(files)
    filename{i} = files(i).name;
    temp = strsplit(filename{i},'_');
    %time{i} = temp{1}; temp = strsplit(temp{2},'_');
    condition{i} = temp{1};
end
%files2 = dir('*0.01_b_trimmedGenotype082416*.mat');

dof = zeros(length(files),1);
dev = zeros(length(files),1);
count = 1;
outputFile = fopen([pathname,'causalVariantList_stringent.txt'],'w');
fprintf(outputFile,['Time\tCondition\tvariantPos(FM)\tvariantPos(FW)\thistory\t',...
    'sum_abs_x(FW)\tdev\tdof\tdev_lenient\tdev_stringent\tchrom\tPVAL\tpos_fineMapping\tpos_fwselection\t',...
    'peakWidth\tpeakProminance\tb_fwselection\tvarianceExplained\tLOD_pos_fineMapping\tLOD_pos_fwselection',...
    '\tref\talt\tRMgenotype\tYJMgenotype',...
    '\tfull_annotation\tclassifier\tmutation\tORF\tORFannotation\n']);

for fileIndex = 1:length(files)
    clear minAlt; clear maxP; clear minPh2; clear fwP; clear resolution; clear pks; %%% These variables are not indexed by fileIndex so they need to be cleared.
    load([pathname,filename{fileIndex}]);
    history = b.history;
    load([pathname,files(fileIndex).name]);
    b.history = history;
    if length(b.history) ~= length(b.bPos)
        fprintf([num2str(length(b.history)),',',num2str(fileIndex),'\n']);
    end
    %%% Don't write the normrnd and modelNoise files
    if ~strcmp(time{fileIndex},'0h')
        
        %%% Calculate percentage of variance explained by each predictor in
        %%% the model
        sumR = zeros(length(b.bPos),1);
        varianceExplained = zeros(length(b.bPos),1);
        for i = 1:length(b.bPos)
            newResidual = b.stats.yr + b.b_fwselection(b.bPos(i))*b.x(:,b.bPos(i));
            sumR(i) = sum(newResidual.^2) - b.stats.SSresid;
        end
        for i = 1:length(b.bPos)
            varianceExplained(i) = sumR(i)/sum(sumR)*b.dev_fwselection;
        end
        
        %%% Use bootstrap p-values to determine a cutoff for fine mapping.
%         bootstrap_pval = [];
%         for i = 1:length(b.bootstrap_b);
%             bootstrap_pval = [bootstrap_pval; -log10(b.bootstrap_pval{i}(b.bootstrap_inmodel{i}))];
%         end
%         pval = -log10(b.stats.PVAL(b.inmodel));
%         increment = 0.1;
%         range = 0:increment:max(pval(~isinf(pval)));
%         FP = zeros(length(range),1);
%         for i = 1:length(range);
%             FP(i) = sum(bootstrap_pval > range(i))/length(b.bootstrap_b)/sum(pval > range(i));
%         end
%         %%%% Cutoff by running a full stepwise selection on nBootstrap
%         %%%% permutations and looking for a threshold where the # of p-val >
%         %%%% threshold in the real model is 20x greater than # of p-val >
%         %%%% threshold in an average bootstrap
%         if ~isempty(range(min(find(FP<0.05))))
%             cutoff(fileIndex) = range(min(find(FP<0.05)));
%         else
%             cutoff(fileIndex) = max(range);
%         end
%         % Use the max bootstrap value from permutation test at each step of the
%         % stepwise selection
%         recursiveBootstrap = zeros(length(b.bootstrap_stepwise),1);
%         for i = 1:length(b.bootstrap_stepwise)
%             recursiveBootstrap(i) = -log10(b.bootstrap_stepwise{i}(round(0.05*length(b.bootstrap_stepwise{i}))));
%         end
%         cutoff2(fileIndex) = max(recursiveBootstrap);
%         nVariants = sum(pval>cutoff(fileIndex));
        
        %these may need to be reordered to reflect order files are loaded...
        %cutoff =[5.15,6.73,6.19,5.31,5.33,5.84,5.02,5];
        %cutoff=[6.73,6.19,5.31,5.15]; %reorder to alphabetical
        cutoff=ones(1,28)*4;
        %cutoff=[10 10 13 14 15 12 8 8 11 10 11 15 10 13 10 10 8 13 9 12 11 15 14 14 8 10 14 10];
        
        bPos_filtered = find(-log10(b.stats.PVAL) > cutoff(fileIndex));
        [~,~,r_filtered] = regress(b.y,[ones(length(b.y),1),b.x(:,bPos_filtered)]);
        
        startIndex = sum(bPos_filtered <= 13) + 1;
        
        dof(fileIndex) = sum(b.bPos > 13);
        dev(fileIndex) = b.dev_fwselection;
        devLenient = sum(varianceExplained(-log10(b.stats.PVAL(b.bPos)) > cutoff(fileIndex)));
        devStringent = sum(varianceExplained(-log10(b.stats.PVAL(b.bPos)) > cutoff(fileIndex)));
        
        if isfield(b,'ph2')
            %         if length(b.bPos) > 200
            %             [~,sortedBpos] = sort(b.bMagLod);
            %             b.bPos = b.bPos(sortedBpos(length(sortedBpos)-199:length(sortedBpos)));
            %             b.bPos = sort(b.bPos);
            %         end
            resolution = zeros(length(bPos_filtered),1);
            pks = zeros(length(bPos_filtered),1);
            fwP = zeros(length(bPos_filtered),1);
            minPh2 = cell(length(bPos_filtered),1);
            fprintf(['fileIndex',num2str(fileIndex),'\n']);
            fprintf(['bPos',num2str(length(bPos_filtered)),'\n']);
            for posIndex = startIndex:length(bPos_filtered);
                pos1 = bPos_filtered(posIndex);
                tempVar = strsplit(variantPos{pos1-13},':');
                chrom = tempVar{1};
                position = str2num(tempVar{2});
                ref = tempVar{3};
                alt = tempVar{4};
                if strcmp(tempVar{5},'0');
                    RMgenotype = 'ref';
                elseif strcmp(tempVar{5},'1');
                    RMgenotype = 'alt';
                else
                    RMgenotype = tempVar{5};
                end
                if strcmp(tempVar{6},'0');
                    YJMgenotype = 'ref';
                elseif strcmp(tempVar{6},'1');
                    YJMgenotype = 'alt';
                else
                    YJMgenotype = tempVar{6};
                end
                % Scan all positions up to 10 kb in either direction for possible
                % causal variants that fit better by the allele swap criteria.
                lower = min(find(strcmp(chrom,allChrom) & allPos > position - 10000)+13);
                upper = max(find(strcmp(chrom,allChrom) & allPos < position + 10000)+13);
                nPosInModel{fileIndex}(posIndex) = sum(b.bPos >= lower & b.bPos <= upper);
                posInModel = b.bPos(find(b.bPos >= lower & b.bPos <= upper));
                queryPos = bPos_filtered(posIndex);
                fwP(posIndex) = -log10(b.stats.PVAL(pos1));
                
                if fwP(posIndex) > cutoff(fileIndex)
                    %%% FilteredPos is a subset of bPos, but b.ph2 is indexed by
                    %%% bPos, so we need a conversion index.
                    filteredIndex = find(pos1==b.bPos);
                    
                    %%% Find the minimum of each row of the Ph2 matrix (Ph2 is
                    %%% the likelihood of the alternate hypothesis H2, and the
                    %%% maximum of the Ph1 (Likelihood of H1).
                    leftBound = zeros(length(b.ph2{posIndex}),1);
                    rightBound = zeros(length(b.ph2{posIndex}),1);
                    minPh2{posIndex} = zeros(length(b.ph2{posIndex}),1);
                    degeneracyCheck = zeros(length(b.ph2{posIndex}),1);
                    for i = 1:length(b.ph2{posIndex});
                        if (sum(b.x(:,lower+i-1) == 1) > length(b.y)/4) && (sum(b.x(:,lower+i-1) == -1) > length(b.y)/4)
                            degeneracyCheck(i) = 0;
                        else
                            degeneracyCheck(i) = 1;
                        end
                    end
                    
                    for i = 1:length(b.ph2{posIndex})
                        temp2 = b.ph2{posIndex}(i,:);
                        temp2(i) = 11; temp2(temp2 == -1) = 11; temp2(temp2 == 0) = 11; %temporary value so as not to affect min calculation, p-value is never 11 or -1;
                        % If we landed in a degenerate row (i+1 and i-1 positions
                        % are -1, then set to zero.
                        if min(degeneracyCheck(max(1,i-1):min(length(temp2),i+1)))==1
                            minPh2{posIndex}(i) = 0;
                            leftBound(i) = 1; rightBound(i) = length(temp2);
                        else
                            [~,leftBound(i)] = min(temp2(1:i));
                            [tempMin,rightBound(i)] = min(temp2(i:length(temp2)));
                            rightBound(i) = i + rightBound(i) - 1;
                            rightBound(i) = max(rightBound(i),floor((i+length(temp2))/2));
                            leftBound(i) = min(leftBound(i),ceil((i+1)/2));
                            
                            if tempMin == 11
                                rightBound(i) = length(temp2);
                            end
                            temp2(i) = min(temp2); temp2(temp2 == 11) = min(temp2); % new temporary value so as not to affect min(-log10(temp2)), p-value is never 11 or -1;
                            if min(temp2) ~= 11
                                minPh2{posIndex}(i) = min(-log10(temp2(leftBound(i):rightBound(i))));
                            else
                                minPh2{posIndex}(i) = 0;
                            end
                        end
                    end
                    
                    queryLoc = pos1 - lower + 1;
                    
                    %%% If the position from forward selection falls on a
                    %%% degenerate column of x that does not give a good
                    %%% fineMappingLod readout (usually lots of missing data or
                    %%% only 1 dominant genotype at the position
                    if max(b.ph2{posIndex}(queryLoc,:)) == -1
                        primaryLoc = queryLoc;
                        resolution(posIndex) = 1;
                        pks(posIndex) = 0;
                    else
                        [pks(posIndex),primaryLoc] = max(minPh2{posIndex}(leftBound(queryLoc):rightBound(queryLoc)));
                        primaryLoc = primaryLoc + leftBound(queryLoc) - 1;
                        fineMappingCutoff = 0.05;
                        % Find the 95% confidence interval for the causal variant
                        % If all other alternative hypotheses have a likelihood <
                        % 0.05, then, we have a single causal variant
                        if pks(posIndex) > -log10(fineMappingCutoff);
                            locs = find(minPh2{posIndex}(leftBound(queryLoc):rightBound(queryLoc)) > pks(posIndex)*0.7);
                            locs = locs + leftBound(queryLoc) - 1;
                            
                            tempVar = strsplit(variantPos{lower+min(locs)-1-13},':');
                            lowerpos = str2num(tempVar{2});
                            tempVar = strsplit(variantPos{lower+max(locs)-1-13},':');
                            upperpos = str2num(tempVar{2});
                            resolution(posIndex) = upperpos-lowerpos+1;
                            % Otherwise, find the interval (set of variants) such that
                            % all other sites have p<0.05 relative to the set.
                        elseif pks(posIndex) > 0
                            temp2 = b.ph2{posIndex}(primaryLoc,:);
                            temp2(temp2==-1) = min(temp2(temp2>0));
                            locs = [primaryLoc,find(-log10(temp2) < -log10(fineMappingCutoff))];
                            locs = unique(locs);
                            temp2(min(locs):max(locs)) = min(temp2(temp2>0));
                            [pks(posIndex)] = min(-log10(temp2(leftBound(primaryLoc):rightBound(primaryLoc))));
                            
                            % Empirical finding is that queryLoc performs better in
                            % model datasets than the location found from pks.
                            primaryLoc = queryLoc;
                            
                            tempVar = strsplit(variantPos{lower+min(locs)-1-13},':');
                            lowerpos = str2num(tempVar{2});
                            tempVar = strsplit(variantPos{lower+max(locs)-1-13},':');
                            upperpos = str2num(tempVar{2});
                            resolution(posIndex) = upperpos-lowerpos+1;
                            
                            %%%% Iterating this process of selection is
                            %%%% pathological because transitive relationships are
                            %%%% broken - eventually all sites are included. n=1 is
                            %%%% good enough, since sites not in the 1st iteraton
                            %%%% are ruled out by exlusion from
                            %%%% b.ph2{posIndex}(locs,:)
                            %                     whileCondition = 1;
                            %                     while whileCondition
                            %                         newLocs = [];
                            %                         for locIndex = 1:length(locs)
                            %                             temp2 = b.ph2{posIndex}(locs(locIndex),:);
                            %                             temp2(temp2==-1) = min(temp2(temp2>0));
                            %                             newLocs = [newLocs,locs(locIndex),find(-log10(temp2) < -log10(fineMappingCutoff))];
                            %                         end
                            %                         newLocs = unique(newLocs);
                            %                         if length(newLocs) == length(locs)
                            %                             whileCondition = 0;
                            %                         else
                            %                             locs = newLocs;
                            %                         end
                            %                     end
                        else %pathological edge case
                            tempVar = strsplit(variantPos{lower-13},':');
                            lowerpos = str2num(tempVar{2});
                            tempVar = strsplit(variantPos{upper-13},':');
                            upperpos = str2num(tempVar{2});
                            resolution(posIndex) = upperpos-lowerpos+1;
                            primaryLoc = queryLoc;
                            fprintf(['pks<0_fileIndex',num2str(fileIndex),'_posIndex',num2str(posIndex),'\n']);
                        end
                    end
                    
                    FineMappingloc = lower + primaryLoc - 1;
                    tempVar = strsplit(variantPos{FineMappingloc-13},':');
                    FineMappingPosition = str2num(tempVar{2});
                    fprintf(outputFile,['%s\t%s\t%d\t%d\t%d\t%d\t%f\t%d\t%f\t%f\t%s\t%f',...
                        '\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%s\t%s\t',...
                        '%s\t%s\t%s\t%s\t%s\n'],...
                        time{fileIndex},condition{fileIndex},FineMappingloc,pos1,b.history(posIndex),sum(abs(b.x(:,pos1))),...
                        dev(fileIndex),dof(fileIndex),devLenient,devStringent,chrom,fwP(posIndex),...
                        FineMappingPosition,position,resolution(posIndex),pks(posIndex),b.b_fwselection(pos1),...
                        varianceExplained(filteredIndex),b.lod1D(FineMappingloc),b.lod1D(pos1),...
                        ref,alt,RMgenotype,YJMgenotype,...
                        annotation{FineMappingloc-13},classifier{FineMappingloc-13},...
                        mutation{FineMappingloc-13},ORF{FineMappingloc-13},ORFannotation{FineMappingloc-13});
                    count = count + 1;
                end
            end
        end
    end
end
fclose(outputFile);




