%trim results based on FDR cutoffs
clear all

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',2)

blueTriplet=[43 172 226]./256;
orangeTriplet=[248 149 33]./256;

figureCounter=1;

load('nutrientFilename.mat')

vFDR=[5 4 5 4 4 5 5 4 4 5 5 5 5 5 5 5 4 5 4 5 4 5 4 5 5 4 4 5];

inputData=readtable('causalVariantList_stringent.txt');

outputTable=[];

for i=1:length(filename)
    
    vIdx1=ismember(inputData.Condition,filename{i});
    vIdx2=inputData.PVAL>vFDR(i);
    
    outputTable=[outputTable;inputData(logical(vIdx1.*vIdx2),:)];
    
end

for i=1:height(outputTable)
    
    tempRow=outputTable(i,:);
    
    outputTable.Condition{i}=tempRow.Condition{1}(1:(end-4));
    outputTable.Time{i}=tempRow.Condition{1}((end-3):end);
    
end

writetable(outputTable,'causalVariantList_trimmed.csv');





