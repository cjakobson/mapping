%scrape data from mapping output and save table of results

clear all

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',2)
figureCounter=1;


%load info on variants
variantInfo=readtable('/Users/cjakobson/Dropbox/JaroszLab/211028_SegregantProteomicsData_V1/variantInfoStructure.csv');


load('radRemapFilename.mat')

for i=1:length(filename)
    tempStr=strsplit(filename{i},{'-','+','_delta'});
    times{i}=tempStr{1}(1:3);
    conditions{i}=tempStr{1}(5:end);
end

conditions=unique(conditions);
times=unique(times);

load('pValCutoffsFromPerm_FDR_0.05.mat')

o=1;
for q=1:length(conditions)
    
    for l=1:length(times)
          
        load('radRemapFilename.mat')
        pValThresh=vCutoff(ismember(filename,[times{l} ' ' conditions{q} '_delta']));
        
        ratioFilename=['radRemap/' times{l} ' ' conditions{q} '_delta.mat'];
        load(ratioFilename)
        
        %filter by pVal
        bPos=bPos(pValues>pValThresh);
        varianceExplained=varianceExplained(pValues>pValThresh);
        pValues=pValues(pValues>pValThresh);

        isQtn=zeros(length(posToMap),1);
        isQtg=zeros(length(posToMap),1);
        hasCandidate=zeros(length(posToMap),1);
        bestCandidate=zeros(length(posToMap),1);
        nCandidates=zeros(length(posToMap),1);
        for m=1:length(isQtn)

            vTemp=candidates{m};
            nCandidates(m)=length(vTemp);

            if length(vTemp)==1
                if vTemp==11    %QTL confirmed as QTN
                    isQtn(m)=1;
                    isQtg(m)=0;
                    hasCandidate(m)=1;
                    bestCandidate(m)=posToMap(m)+vTemp-11;
                else            %QTN shifted relative to QTL
                    hasCandidate(m)=1;
                    bestCandidate(m)=posToMap(m)+vTemp-11;
                    isQtn(m)=1;
                    isQtg(m)=0;
                end
            else %check for QTGs
                tempIdx=posToMap(m)+vTemp-11;
                if (length(tempIdx)>1)%&&(sum(ismember(vTemp,11))>0)
                    toLookup=mod(tempIdx-41,12054);
                    tempTable=array2table([variantInfo.gene1(toLookup)...
                        variantInfo.gene2(toLookup)]);
                    tempUnique=unique(tempTable,'rows');
                    if height(tempUnique)==1
                        isQtn(m)=0;
                        isQtg(m)=1;
                        hasCandidate(m)=0;
                        bestCandidate(m)=0;
                    end
                end
            end

        end

        qtnPos=posToMap(logical(isQtn));
        qtgPos=posToMap(logical(isQtg));
        candidatePos=posToMap(logical(hasCandidate));

        qtlIsQtn=ismember(bPos,qtnPos);
        qtlIsQtg=ismember(bPos,qtgPos);
        qtlHasCandidate=ismember(bPos,candidatePos);

        for m=1:length(bPos)

            vCondition{o}=conditions{q};
            vTime{o}=times{l};

            vBpos(o)=bPos(m);
            vPval(o)=pValues(m);
            vBeta(o)=b_fwselection(bPos(m));
            vVarExp(o)=varianceExplained(m);

%             if bPos(m)>41
%                 vIndex(o)=mod(bPos(m)-41,12054);
%             else
%                 vIndex(o)=0;
%             end

            if qtlIsQtn(m)
                vIsQtn(o)=1;
            else
                vIsQtn(o)=0;
            end
            
            if qtlIsQtg(m)
                vIsQtg(o)=1;
            else
                vIsQtg(o)=0;
            end
            
            if qtlHasCandidate(m)
                vFineCandidate(o)=bestCandidate(ismember(posToMap,bPos(m)));
            else
                vFineCandidate(o)=0;
            end
            
            if sum(ismember(posToMap,bPos(m)))>0
                vNcandidates(o)=nCandidates(ismember(posToMap,bPos(m)));
            end

            o=o+1;

        end

    end
    
end



vFineCandidate(vFineCandidate==0)=vBpos(vFineCandidate==0);

vIndex(vFineCandidate>41)=mod(vFineCandidate(vFineCandidate>41)-41,12054);
vIndex(vFineCandidate<=41)=0;

toOutput=table(vCondition',vTime',vBpos',vPval',vBeta',vVarExp',vIsQtn',vIsQtg',vFineCandidate',vNcandidates',...
    'VariableNames',{'condition','time','bPos','pVal','beta','varExp',...
    'isQtn','isQtg','bestCandidate','nCandidates'});

%make dummy info row for geometric terms
variantInfo{height(variantInfo)+1,7}={'geometric'};

vIndex(vIndex==0)=height(variantInfo);



toOutput=[toOutput variantInfo(vIndex,:)];



%load genotype info
load('radRemapFilename.mat')
load('radRemapTrait.mat')

load('phasedVLCgenotype.mat')
genotypes=phasedVLCgenotype;
clear phasedVLCgenotype

nPlates=length(trait{1})/384;
nStrains=length(trait{1});

genotypes=genotypes(1:nStrains,:);

[~,temp]=size(genotypes);
nLoci=temp/4;    

modelGenotypes=zeros(nStrains,2*nLoci);
for i=1:nStrains
    vGenotype=[genotypes(i,1:nLoci)-genotypes(i,(nLoci+1):(2*nLoci)) genotypes(i,(2*nLoci+1):(3*nLoci))+genotypes(i,(3*nLoci+1):(4*nLoci))];
    modelGenotypes(i,:)=vGenotype;
end

clear genotypes;


m=1;
%calculate delta Z-score to call buffer/potentiated/inverted
for i=1:length(conditions)
    
    for j=1:length(times)
    
        baselineQuery=[times{j} ' ' conditions{i} '-rad'];
        baselineTrait=trait{ismember(filename,baselineQuery)};

        bufferQuery=[times{j} ' ' conditions{i} '+rad'];
            
        if sum(ismember(filename,bufferQuery))>0

            bufferTrait=trait{ismember(filename,bufferQuery)};

            conditionIdx=ismember(toOutput.condition,conditions{i});
            timeIdx=ismember(toOutput.time,times{j});

            tempIdx=logical(conditionIdx.*timeIdx);

            %tempBpos=toOutput.bPos(tempIdx);
            tempBpos=toOutput.bestCandidate(tempIdx);

            for l=1:length(tempBpos)

                posQuery=tempBpos(l)-41;

                if (posQuery<=12054)&&(posQuery>0)

                    idx1=modelGenotypes(:,posQuery)==1;
                    idx2=modelGenotypes(:,posQuery)==-1;

                    deltaZbaseline(m)=mean(baselineTrait(idx1),'omitnan')-...
                        mean(baselineTrait(idx2),'omitnan');
                    deltaZbuffer(m)=mean(bufferTrait(idx1),'omitnan')-...
                        mean(bufferTrait(idx2),'omitnan');
                    
                    %dominanceBaseline(m)=0;
                    %dominanceBuffer(m)=0;
                    
                    %calculate dominance for linear loci too, for
                    %prioritization
                    %also calcultate extent of dominance [fully RM = 1;
                    %fully YJM = -1]
                    %idx3=modelGenotypes(:,posQuery-12054)==1;
                    idx3=modelGenotypes(:,posQuery+12054)==0;
                    %idx4=modelGenotypes(:,posQuery-12054)==-1;
                    idx4=modelGenotypes(:,posQuery+12054)==1;
                    
                    tempRmZ=mean(baselineTrait(idx3),'omitnan');
                    tempYjmZ=mean(baselineTrait(idx4),'omitnan');
                    
                    %quantities are swapped relative to below
                    if (deltaZbaseline(m)*tempRmZ)>0
                        %dominanceBaseline(m)=deltaZbaseline(m)/tempRmZ;
                        dominanceBaseline(m)=tempRmZ/deltaZbaseline(m);
                    elseif (deltaZbaseline(m)*tempYjmZ)>0
                        %dominanceBaseline(m)=-deltaZbaseline(m)/tempYjmZ;
                        dominanceBaseline(m)=-tempYjmZ/deltaZbaseline(m);
                    else
                        dominanceBaseline(m)=0;
                    end
                    
                    tempRmZ=mean(bufferTrait(idx3),'omitnan');
                    tempYjmZ=mean(bufferTrait(idx4),'omitnan');
                    
                    if (deltaZbuffer(m)*tempRmZ)>0
                        %dominanceBuffer(m)=deltaZbuffer(m)/tempRmZ;
                        dominanceBuffer(m)=tempRmZ/deltaZbuffer(m);
                    elseif (deltaZbuffer(m)*tempYjmZ)>0
                        %dominanceBuffer(m)=-deltaZbuffer(m)/tempYjmZ;
                        dominanceBuffer(m)=-tempYjmZ/deltaZbuffer(m);
                    else
                        dominanceBuffer(m)=0;
                    end

                elseif posQuery>12054

                    idx1=modelGenotypes(:,posQuery)==0;
                    idx2=modelGenotypes(:,posQuery)==1;

                    deltaZbaseline(m)=mean(baselineTrait(idx1),'omitnan')-...
                        mean(baselineTrait(idx2),'omitnan');
                    deltaZbuffer(m)=mean(bufferTrait(idx1),'omitnan')-...
                        mean(bufferTrait(idx2),'omitnan');
                    
                    %also calcultate extent of dominance [fully RM = 1;
                    %fully YJM = -1]
                    idx3=modelGenotypes(:,posQuery-12054)==1;
                    idx4=modelGenotypes(:,posQuery-12054)==-1;
                    
                    tempRmZ=mean(baselineTrait(idx3),'omitnan');
                    tempYjmZ=mean(baselineTrait(idx4),'omitnan');
                    
                    if (deltaZbaseline(m)*tempRmZ)>0
                        dominanceBaseline(m)=deltaZbaseline(m)/tempRmZ;
                    elseif (deltaZbaseline(m)*tempYjmZ)>0
                        dominanceBaseline(m)=-deltaZbaseline(m)/tempYjmZ;
                    else
                        dominanceBaseline(m)=0;
                    end
                    
                    tempRmZ=mean(bufferTrait(idx3),'omitnan');
                    tempYjmZ=mean(bufferTrait(idx4),'omitnan');
                    
                    if (deltaZbuffer(m)*tempRmZ)>0
                        dominanceBuffer(m)=deltaZbuffer(m)/tempRmZ;
                    elseif (deltaZbuffer(m)*tempYjmZ)>0
                        dominanceBuffer(m)=-deltaZbuffer(m)/tempYjmZ;
                    else
                        dominanceBuffer(m)=0;
                    end
                    
                
                end

                m=m+1;

            end

        end
        
    end
    
end


bufferType=cell(length(deltaZbaseline),1);
for i=1:length(deltaZbaseline)
    
    delta1=deltaZbaseline(i);
    delta2=deltaZbuffer(i);
    
    if (delta1*delta2)<0
        
        bufferType{i}='line_crossing';
        
    elseif abs(delta1)<abs(delta2)
        
        bufferType{i}='buffered';
        
    elseif abs(delta1)>abs(delta2)
        
        bufferType{i}='potentiated';
        
    end
    
end


deltaZtable=table(deltaZbaseline',deltaZbuffer',bufferType,dominanceBaseline',dominanceBuffer',...
    'VariableNames',...
    {'deltaZbaseline','deltaZbuffer','interactionType','dominanceBaseline','dominanceBuffer'});

toOutput=[toOutput deltaZtable];


%convert deltaZ into %change in growth using mean and std info
rawGrowthData=readtable('radMeanStdData.csv');

percentBaseline=nan(height(toOutput),1);
percentBuffer=percentBaseline;
for i=1:height(toOutput)
    
    tempCondition=toOutput.condition{i};
    tempTime=toOutput.time{i};
    
    tempQuery=[tempTime ' ' tempCondition '-rad'];
    
    tempIdx=find(ismember(rawGrowthData.conditions,tempQuery));
    
    if sum(tempIdx)>0
        
        mean1=rawGrowthData.Var2(tempIdx);
        mean2=rawGrowthData.Var2(tempIdx+1);
        
        std1=rawGrowthData.Var3(tempIdx);
        std2=rawGrowthData.Var3(tempIdx+1);
        
        percentBaseline(i)=(toOutput.deltaZbaseline(i)*std1)/mean1;
        percentBuffer(i)=(toOutput.deltaZbuffer(i)*std2)/mean2;
        
    end
    
end

percentTable=table(percentBaseline,percentBuffer);

toOutput=[toOutput percentTable];







%add whether Gong client, SGD genetic or physical interactor, etc
tempInput=readtable('GongTableS2.xlsx');

chapLabels={'Hsp82','Hsc82','Hsp104','Hsp42','Ssa1','Ssa2','Ssa3','Ssa4',...
    'Sse1','Ydj1'};

clear upper
for i=1:length(chapLabels)
    
    tempData=table2array(tempInput(:,ismember(tempInput.Properties.VariableNames,chapLabels{i})));
    
    for j=1:length(tempData)
        
        if length(tempData{j})>0
            
            tempStr=strsplit(tempData{j},',');
            
            %make upper case to match bufferInput
            clientArray{i,j}=upper(tempStr{1});
            
        end
        
    end
    
end


for i=1:height(toOutput)
    
    if ~isempty(toOutput.common1{i})
        tempStr=strsplit(toOutput.common1{i},';');
        tempCommon1{i}=tempStr{1};
    end
    if ~isempty(toOutput.common2{i})
        tempStr=strsplit(toOutput.common2{i},';');
        tempCommon2{i}=tempStr{1};
    end
    
end

tempCommon1(cellfun(@isempty,tempCommon1))={'NA'};
tempCommon2(cellfun(@isempty,tempCommon2))={'NA'};

for j=1:length(chapLabels)
        
    tempArray=clientArray(j,:);
    tempArray(cellfun(@isempty,tempArray))=[];

    v1=ismember(tempCommon1,tempArray);
    v2=ismember(tempCommon2,tempArray);

    vChap(:,j)=logical(v1+v2);
        
end


chapTable=array2table(vChap,'VariableNames',chapLabels);

toOutput=[toOutput chapTable];


%also annotate with ase information
aseData=readtable('/Users/cjakobson/Dropbox/JaroszLab/hsp90mapping/harmonizeRnaSeqAnalysis/radAseData.csv');



for i=1:height(toOutput)
    
    if toOutput.index(i)>0

        toOutput2(i,:)=[toOutput(i,:) aseData(toOutput.index(i),[3:6 36:51])];
        
    else
        
        toOutput2(i,:)=[toOutput(i,:) array2table(zeros(1,20))];
        
    end
    
end







writetable(toOutput2,'linearDelta.csv')






