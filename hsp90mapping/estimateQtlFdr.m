%analyze QTL-level FDR from permutations
clear

figureCounter=1; %figure counter

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',2)


load('radRemapFilename.mat')

for n=1:length(filename)
   
    load('radRemapFilename.mat')
    toGet=['radRemap/' filename{n} '.mat'];
    if exist(toGet)
        load(toGet)
        trueP{n}=pValues;
    end
    
    load('radRemapFilename.mat')
    toGet=['radRemap/' filename{n} '_perm.mat'];
    if exist(toGet)
        load(toGet)
        for o=1:length(pValues)
            permP{n,o}=pValues{o};
        end
    end
    
end


[nProteins,nPerms]=size(permP);
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1)
hold on
pValRange=3:0.1:20;
desiredFDR=0.05;
for i=1:nProteins
    
    for j=1:length(pValRange)
    
        tempThresh=pValRange(j);
        
        trueSum(i,j)=sum(trueP{i}>tempThresh);
        
        permSum(i,j)=0;
        for k=1:nPerms
            
            permSum(i,j)=permSum(i,j)+sum(permP{i,k}>tempThresh);
            
        end
    
    end
    
    fdrMat(i,:)=(permSum(i,:)/nPerms)./trueSum(i,:);
    fdrMat(isnan(fdrMat))=Inf;
    
    scatter(pValRange,fdrMat(i,:),10,'k','filled')
    
    %calculate threshold to constrain FDR
    tempIdx=min(find(fdrMat(i,:)<desiredFDR));
    if isempty(tempIdx)
        tempIdx=length(pValRange);
    end
    
    vCutoff(i)=pValRange(tempIdx);
    vFdr(i)=fdrMat(i,tempIdx);
    
end


plot(pValRange,median(fdrMat,'omitnan'),'r')
ylim([0 1])


subplot(2,3,2)
histogram(vFdr)
xlim([0 0.1])

subplot(2,3,3)
histogram(vCutoff)


save(['pValCutoffsFromPerm_FDR_' num2str(desiredFDR) '.mat'],'vCutoff')


