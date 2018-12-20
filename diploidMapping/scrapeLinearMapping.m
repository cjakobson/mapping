%scrape data from fitting output and save table of results

clear all

%load info on variants
inputMat=tdfread('variantPos080316_trimmedGenotype.txt');

%load TSS info to assign intergenic variants to nearest gene TSS
%load Snyder RNAseq TSS data
[num,txt,~]=xlsread('TableS4.xls');
tssGenes=txt(4:end,1);
tempTssChr=txt(4:end,2);
tssPos=num(:,3);
tesPos=num(:,4);
tssIsWatson=num(:,1)<num(:,2);

chromArray={'I','II','III','IV','V','VI','VII','VIII','IX','X',...
    'XI','XII','XIII','XIV','XV','XVI','Mito','chrMito'};

%recast chromsome as num
for i=1:length(tempTssChr)
    
    tssChr(i)=find(ismember(chromArray,tempTssChr{i}(4:end)));
    
end


%load S288C (ref) sequence
refSeq=fastaread('S288C_reference_sequence_R64-2-1_20150113.fsa');

chromArray2={'I      ','II     ','III    ','IV     ','V      ','VI     ',...
    'VII    ','VIII   ','IX     ','X      ',...
    'XI     ','XII    ','XIII   ','XIV    ',...
    'XV     ','XVI    ','chrMito'};

%read out chromosome and positions for all extragenic variants in the
%cross as well as genes and positions for all syn variants

%first categorize all segregating variants
allVariantGenes=cell(length(inputMat.INFO),1);
allVariantType=cell(length(inputMat.INFO),1);
synVariantPos=cell(length(inputMat.INFO),1);

for i=1:length(inputMat.INFO)
    
    variantInfo=inputMat.INFO(i,:);
                    
    strArray=strsplit(variantInfo,'|');
    allVariantGenes{i}=strArray{4};
    allVariantType{i}=strArray{2};
    if length(strArray)>9
    	synVariantPos{i}=strArray{10};
    end

end

%also need ref and alt codons for all segregating variants...
%find ref and alt codons

synonymousIdx=ismember(allVariantType,'synonymous_variant');

synLookup=find(synonymousIdx);
for i=1:length(synLookup)
    variantIdx=synLookup(i);
    variantInfo=inputMat.INFO(variantIdx,:);
                    
    strArray=strsplit(variantInfo,'|');
    tempPos=strArray{10};
                    
    tempGenomePos=inputMat.POS(variantIdx,:);
    tempChr=inputMat.CHROM(variantIdx,:);
    tempOrfPos=str2num(tempPos(3:(end-3)));
    codonPos=mod(tempOrfPos,3);

    %check for Watson/Crick to grab correct bases
    tempGeneName=strArray{4};
    
    if length(tempGeneName)>6
        tempWatsonCrick=tempGeneName(7);
    end
    
    refBase=strArray{10}(end-2);
    altBase=strArray{10}(end);

    if tempWatsonCrick=='W'
        if codonPos==0
            allRefCodon{i}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-2):(tempGenomePos));
            allAltCodon{i}=allRefCodon{i};
            allRefCodon{i}(3)=refBase;
            allAltCodon{i}(3)=altBase;
            
            allRefDicodon1{i}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-5):(tempGenomePos));
            allRefDicodon2{i}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-2):(tempGenomePos+3));
        elseif codonPos==1
            allRefCodon{i}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos):(tempGenomePos+2));
            allAltCodon{i}=allRefCodon{i};
            allRefCodon{i}(1)=refBase;
            allAltCodon{i}(1)=altBase;
            
            allRefDicodon1{i}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-3):(tempGenomePos+2));
            allRefDicodon2{i}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos):(tempGenomePos+5));
        elseif codonPos==2
            allRefCodon{i}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-1):(tempGenomePos+1));
            allAltCodon{i}=allRefCodon{i};
            allRefCodon{i}(2)=refBase;
            allAltCodon{i}(2)=altBase;
            
            allRefDicodon1{i}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-4):(tempGenomePos+1));
            allRefDicodon2{i}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos+1):(tempGenomePos+4));
        else
            allRefCodon{i}=[];
            allAltCodon{i}=[];
            
            allRefDicodon1{i}=[];
            allRefDicodon2{i}=[];
        end
    elseif tempWatsonCrick=='C'
        if codonPos==0
            allRefCodon{i}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos):(tempGenomePos+2)));
            allAltCodon{i}=allRefCodon{i};
            allRefCodon{i}(3)=refBase;
            allAltCodon{i}(3)=altBase;
            
            allRefDicodon1{i}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos):(tempGenomePos+5)));
            allRefDicodon2{i}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-3):(tempGenomePos+2)));
        elseif codonPos==1
            allRefCodon{i}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-2):(tempGenomePos)));
            allAltCodon{i}=allRefCodon{i};
            allRefCodon{i}(1)=refBase;
            allAltCodon{i}(1)=altBase;
            
            allRefDicodon1{i}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-5):(tempGenomePos)));
            allRefDicodon2{i}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-2):(tempGenomePos+3)));
        elseif codonPos==2
            allRefCodon{i}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-1):(tempGenomePos+1)));
            allAltCodon{i}=allRefCodon{i};
            allRefCodon{i}(2)=refBase;
            allAltCodon{i}(2)=altBase;
            
            allRefDicodon1{i}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-1):(tempGenomePos+4)));
            allRefDicodon2{i}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-4):(tempGenomePos+1)));
        else
            allRefCodon{i}=[];
            allAltCodon{i}=[];
            
            allRefDicodon1{i}=[];
            allRefDicodon2{i}=[];
        end
    else
        refCodon{i}=[];
        altCodon{i}=[];
            
        allRefDicodon1{i}=[];
        allRefDicodon2{i}=[];
    end
end

toOutput=table(allRefCodon',allAltCodon','VariableNames',...
{'ref','alt'});
writetable(toOutput,'controlCodons.csv');

toOutput=table(allRefDicodon1',allRefDicodon2','VariableNames',...
{'dicodon1','dicodon2'});
writetable(toOutput,'controlDicodons.csv');



extragenicIdx=logical(ismember(allVariantType,'upstream_gene_variant')+...
    ismember(allVariantType,'downstream_gene_variant')+...
    ismember(allVariantType,'intergenic_region'));

%output chromosomes and positions for extragenic variants
extragenicChrom=inputMat.CHROM(extragenicIdx,:);
extragenicPos=inputMat.POS(extragenicIdx,:);

for i=1:length(extragenicPos)
    
    %find nearest genes upstream and downstream and
    %call promoter status by direction
    tempChr=find(ismember(chromArray2,extragenicChrom(i,:)));

    tempIdx=(tssChr==tempChr);
    tempDist=tssPos(tempIdx)-extragenicPos(i);
    tempGenes=tssGenes(tempIdx);

    toKeep=~isnan(tempDist);
    tempDist=tempDist(toKeep);
    tempGenes=tempGenes(toKeep);

    idxUp=find(tempDist<0);
    idxDown=find(tempDist>0);
    
    if ~isempty(idxUp)&&~isempty(idxDown)&&(tempChr<17) %mito doesn't have W/C
        upGene=tempGenes{idxUp(end)}(7);
        downGene=tempGenes{idxDown(1)}(7);

        if upGene=='W'
            if downGene=='W'
                extragenicPromoter{i}='unidirectional';
            elseif downGene=='C'
                extragenicPromoter{i}='head_on';
            else
                extragenicPromoter{i}=[];
            end
        elseif upGene=='C'
            if downGene=='W'
                extragenicPromoter{i}='divergent';
            elseif downGene=='C'
                extragenicPromoter{i}='unidirectional';
            else
                extragenicPromoter{i}=[];
            end
        else
            extragenicPromoter{i}=[];
        end
    else %end of chr
        extragenicPromoter{i}='unidirectional';
    end
    
end

toOutput=table(extragenicChrom,extragenicPos,extragenicPromoter','VariableNames',...
{'chrom','pos','promoter'});
writetable(toOutput,'extragenicPositions.csv');


%output genes and positions for synonymous variants
synonymousGenes=allVariantGenes(synonymousIdx);
synonymousPos=synVariantPos(synonymousIdx);

toOutput=table(synonymousGenes,synonymousPos,'VariableNames',...
{'gene','pos'});
writetable(toOutput,'synonymousPositions.csv');


%output genes and positions for missense variants
missenseIdx=find(ismember(allVariantType,'missense_variant'));


missenseGenes=allVariantGenes(missenseIdx);

for i=1:length(missenseIdx)
    
    tempIdx=missenseIdx(i);
    variantInfo=inputMat.INFO(tempIdx,:);
                    
    strArray=strsplit(variantInfo,'|');
    
    if length(strArray)>10
        
        missensePosInOrf(i)=str2num(strArray{11}(6:(end-3)));
        
    end
end

toOutput=table(missenseGenes,missensePosInOrf','VariableNames',...
{'gene','posInOrf'});
writetable(toOutput,'missensePositions.csv');


%output positions of all segregating variants
allChrom=inputMat.CHROM;
allPos=inputMat.POS;

toOutput=table(allChrom,allPos,'VariableNames',...
{'chrom','pos'});
writetable(toOutput,'allPositions.csv');



%now address CVs

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',2)

pCutoffQTL=5;
pCutoffQTN=5;

load('filenameMerged.mat')


nConditions=15;

conditions={'2%EtOH','2%Gal','2%Glc 37C','2.7uM clotrimazole','2%Glc','2%Gly','2%Raff','2%Mal',...
    '100mM Ca','2mM Co','100mM Li','100mM Mg','5mM Mn','1mM Ni','2.5mM Zn'};
times={'24h','48h','72h','96h','120h','144h'};

%scrape mapping results for each condition and time point    
m=1;
for n=1:length(filename)
    
    load([filename{n} '.mat'])
    load('filenameMerged.mat')
    
    for l=1:length(bPos)

        %test for pValue
        if pValues(l)>pCutoffQTL
            conditionVector{m}=conditions(ceil(n/6));
            
            if mod(n,6)==0
                timeVector{m}=times(6);
            else
                timeVector{m}=times(mod(n,6));
            end

            bPosVector(m)=bPos(l);
            pValVector(m)=pValues(l);
            
            betaVector(m)=b_fwselection(bPos(l));
            varExpVector(m)=varianceExplained(l);
            
            %test for QTN or not
            if sum(ismember(posToMap,bPos(l)))==1
                
                mappingIdx=find(posToMap==bPos(l));
            
                if (length(candidates{mappingIdx})==1)&&(pValVector(l)>pCutoffQTN)
                    
                    isQTNvector(m)=1;
                    isQTGvector(m)=0;
                    qtnPosVector(m)=posToMap(mappingIdx);
                    
                    %if QTN, also look up gene and variant type
                    variantIdx=mod(qtnPosVector(m)-53,12054);
                    
                    variantInfo=inputMat.INFO(variantIdx,:);
                    
                    strArray=strsplit(variantInfo,'|');
                    geneVector{m}=strArray{4};
                    varTypeVector{m}=strArray{2};
                    
                    if strcmp(varTypeVector{m},'synonymous_variant')
        
                        position{m}=strArray{10};
                        encoded{m}=[];     
                        chr{m}=inputMat.CHROM(variantIdx,:);
                        nucleotide(m)=inputMat.POS(variantIdx,:);
                        
                        %find ref and alt codons
                        tempGenomePos=inputMat.POS(variantIdx,:);
                        tempChr=inputMat.CHROM(variantIdx,:);
                        tempOrfPos=str2num(position{m}(3:(end-3)));
                        codonPos=mod(tempOrfPos,3);
                        %check for Watson/Crick to grab correct bases
                        tempGeneName=strArray{4};

                        if length(tempGeneName)>6
                            tempWatsonCrick=tempGeneName(7);
                        end
                        
                        refBase=strArray{10}(end-2);
                        altBase=strArray{10}(end);
                        
                        if tempWatsonCrick=='W'
                            if codonPos==0
                                refCodon{m}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-2):(tempGenomePos));
                                altCodon{m}=refCodon{m};
                                refCodon{m}(3)=refBase;
                                altCodon{m}(3)=altBase;
                                
                                refDicodon1{m}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-5):(tempGenomePos));
                                refDicodon2{m}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-2):(tempGenomePos+3));
                            elseif codonPos==1
                                refCodon{m}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos):(tempGenomePos+2));
                                altCodon{m}=refCodon{m};
                                refCodon{m}(1)=refBase;
                                altCodon{m}(1)=altBase;
                                
                                refDicodon1{m}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-3):(tempGenomePos+2));
                                refDicodon2{m}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos):(tempGenomePos+5));
                            elseif codonPos==2
                                refCodon{m}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-1):(tempGenomePos+1));
                                altCodon{m}=refCodon{m};
                                refCodon{m}(2)=refBase;
                                altCodon{m}(2)=altBase;
                                
                                refDicodon1{m}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-4):(tempGenomePos+1));
                                refDicodon2{m}=refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos+1):(tempGenomePos+4));
                            else
                                refCodon{m}=[];
                                altCodon{m}=[];
                                
                                refDicodon1{m}=[];
                                refDicodon2{m}=[];
                            end
                        elseif tempWatsonCrick=='C'
                            if codonPos==0
                                refCodon{m}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos):(tempGenomePos+2)));
                                altCodon{m}=refCodon{m};
                                refCodon{m}(3)=refBase;
                                altCodon{m}(3)=altBase;
                                
                                refDicodon1{m}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos):(tempGenomePos+5)));
                                refDicodon2{m}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-3):(tempGenomePos+2)));
                            elseif codonPos==1
                                refCodon{m}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-2):(tempGenomePos)));
                                altCodon{m}=refCodon{m};
                                refCodon{m}(1)=refBase;
                                altCodon{m}(1)=altBase;
            
                                refDicodon1{m}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-5):(tempGenomePos)));
                                refDicodon2{m}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-2):(tempGenomePos+3)));
                            elseif codonPos==2
                                refCodon{m}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-1):(tempGenomePos+1)));
                                altCodon{m}=refCodon{m};
                                refCodon{m}(2)=refBase;
                                altCodon{m}(2)=altBase;
                                
                                refDicodon1{m}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-1):(tempGenomePos+4)));
                                refDicodon2{m}=seqrcomplement(refSeq(ismember(chromArray2,tempChr)).Sequence((tempGenomePos-4):(tempGenomePos+1)));
                            else
                                refCodon{m}=[];
                                altCodon{m}=[];
                                
                                refDicodon1{m}=[];
                                refDicodon2{m}=[];
                            end
                        else
                            refCodon{m}=[];
                            altCodon{m}=[];
                                
                            refDicodon1{m}=[];
                            refDicodon2{m}=[];
                        end
                        
                        promoter{m}=[];
                            
                    elseif strcmp(varTypeVector{m},'missense_variant')

                        position{m}=strArray{10};
                        encoded{m}=strArray{11};   
                        chr{m}=[];
                        nucleotide(m)=0;
                        refCodon{m}=[];
                        altCodon{m}=[];    
                        refDicodon1{m}=[];
                        refDicodon2{m}=[];
                        promoter{m}=[];
                        
                    elseif strcmp(varTypeVector{m},'upstream_gene_variant')
                        
                        varTypeVector{m}='intergenic_region';
                        chr{m}=inputMat.CHROM(variantIdx,:);
                        nucleotide(m)=inputMat.POS(variantIdx,:);
                        position{m}=refSeq(ismember(chromArray2,chr{m})).Sequence((nucleotide(m)-5):(nucleotide(m)+5));
                        encoded{m}=strArray{18};
                        refCodon{m}=[];
                        altCodon{m}=[];
                        refDicodon1{m}=[];
                        refDicodon2{m}=[];
                        
                        %find nearest genes upstream and downstream and
                        %call promoter status by direction
                        tempArray=strsplit(chr{m},' ');
                        tempChr=find(ismember(chromArray,tempArray{1}));
                        
                        tempIdx=(tssChr==tempChr);
                        tempDist=tssPos(tempIdx)-nucleotide(m);
                        tempGenes=tssGenes(tempIdx);
                        
                        toKeep=~isnan(tempDist);
                        tempDist=tempDist(toKeep);
                        tempGenes=tempGenes(toKeep);
                        
                        idxUp=find(tempDist<0);
                        idxDown=find(tempDist>0);
                        
                        if ~isempty(idxUp)&&~isempty(idxDown)
                            upGene=tempGenes{idxUp(end)}(7);
                            downGene=tempGenes{idxDown(1)}(7);

                            if upGene=='W'
                                if downGene=='W'
                                    promoter{m}='unidirectional';
                                elseif downGene=='C'
                                    promoter{m}='head_on';
                                else
                                    promoter{m}=[];
                                end
                            elseif upGene=='C'
                                if downGene=='W'
                                    promoter{m}='divergent';
                                elseif downGene=='C'
                                    promoter{m}='unidirectional';
                                else
                                    promoter{m}=[];
                                end
                            else
                                promoter{m}=[];
                            end
                        else %end of chr
                            promoter{m}='unidirectional';
                        end
                    
                    elseif strcmp(varTypeVector{m},'downstream_gene_variant')

                        varTypeVector{m}='intergenic_region';
                        chr{m}=inputMat.CHROM(variantIdx,:);
                        nucleotide(m)=inputMat.POS(variantIdx,:);
                        position{m}=refSeq(ismember(chromArray2,chr{m})).Sequence((nucleotide(m)-5):(nucleotide(m)+5));
                        encoded{m}=strArray{18};
                        refCodon{m}=[];
                        altCodon{m}=[];
                        refDicodon1{m}=[];
                        refDicodon2{m}=[];
                        
                        %find nearest genes upstream and downstream and
                        %call promoter status by direction
                        tempArray=strsplit(chr{m},' ');
                        tempChr=find(ismember(chromArray,tempArray{1}));
                        
                        tempIdx=(tssChr==tempChr);
                        tempDist=tssPos(tempIdx)-nucleotide(m);
                        tempGenes=tssGenes(tempIdx);
                        
                        toKeep=~isnan(tempDist);
                        tempDist=tempDist(toKeep);
                        tempGenes=tempGenes(toKeep);
                        
                        idxUp=find(tempDist<0);
                        idxDown=find(tempDist>0);
                        
                        if ~isempty(idxUp)&&~isempty(idxDown)
                            upGene=tempGenes{idxUp(end)}(7);
                            downGene=tempGenes{idxDown(1)}(7);

                            if upGene=='W'
                                if downGene=='W'
                                    promoter{m}='unidirectional';
                                elseif downGene=='C'
                                    promoter{m}='head_on';
                                else
                                    promoter{m}=[];
                                end
                            elseif upGene=='C'
                                if downGene=='W'
                                    promoter{m}='divergent';
                                elseif downGene=='C'
                                    promoter{m}='unidirectional';
                                else
                                    promoter{m}=[];
                                end
                            else
                                promoter{m}=[];
                            end
                        else %end of chr
                            promoter{m}='unidirectional';
                        end
                        
                    elseif strcmp(varTypeVector{m},'intergenic_region')
                        
                        chr{m}=inputMat.CHROM(variantIdx,:);
                        nucleotide(m)=inputMat.POS(variantIdx,:);
                        position{m}=refSeq(ismember(chromArray2,chr{m})).Sequence((nucleotide(m)-5):(nucleotide(m)+5));
                        if length(strArray)>=8
                            encoded{m}=strArray{8};
                        else
                            encoded{m}=[];
                        end
                        refCodon{m}=[];
                        altCodon{m}=[];
                        refDicodon1{m}=[];
                        refDicodon2{m}=[];
                        
                        %assign to gene of nearest TSS
                        %convert chr to num
                        tempArray=strsplit(chr{m},' ');
                        tempChr=find(ismember(chromArray,tempArray{1}));
                        
                        %find nearest TSS
                        tempIdx=(tssChr==tempChr);
                        [~,geneIdx]=min(abs(tssPos(tempIdx)-nucleotide(m)));
                        tempGenes=tssGenes(tempIdx);
                        geneVector{m}=tempGenes{geneIdx};
                        
                        %find nearest genes upstream and downstream and
                        %call promoter status by direction
                        tempArray=strsplit(chr{m},' ');
                        tempChr=find(ismember(chromArray,tempArray{1}));
                        
                        tempIdx=(tssChr==tempChr);
                        tempDist=tssPos(tempIdx)-nucleotide(m);
                        tempGenes=tssGenes(tempIdx);
                        
                        toKeep=~isnan(tempDist);
                        tempDist=tempDist(toKeep);
                        tempGenes=tempGenes(toKeep);
                        
                        idxUp=find(tempDist<0);
                        idxDown=find(tempDist>0);
                        
                        if ~isempty(idxUp)&&~isempty(idxDown)
                            upGene=tempGenes{idxUp(end)}(7);
                            downGene=tempGenes{idxDown(1)}(7);

                            if upGene=='W'
                                if downGene=='W'
                                    promoter{m}='unidirectional';
                                elseif downGene=='C'
                                    promoter{m}='head_on';
                                else
                                    promoter{m}=[];
                                end
                            elseif upGene=='C'
                                if downGene=='W'
                                    promoter{m}='divergent';
                                elseif downGene=='C'
                                    promoter{m}='unidirectional';
                                else
                                    promoter{m}=[];
                                end
                            else
                                promoter{m}=[];
                            end
                        else %end of chr
                            promoter{m}='unidirectional';
                        end
                        
                    else
                        
                        position{m}=[];
                        encoded{m}=[];   
                        chr{m}=[];
                        nucleotide(m)=0;
                        refCodon{m}=[];
                        altCodon{m}=[];
                        refDicodon1{m}=[];
                        refDicodon2{m}=[];
                        promoter{m}=[];
                        
                    end
                    
                elseif length(candidates{mappingIdx})>1
                    
                    qtgPosVector=posToMap(mappingIdx)+candidates{mappingIdx}-11;
                    
                    clear qtgGeneVector
                    
                    for x=1:length(qtgPosVector)
                        
                        variantIdx=mod(qtgPosVector(x)-53,12054);
                    
                        variantInfo=inputMat.INFO(variantIdx,:);
                        
                        strArray=strsplit(variantInfo,'|');
                        qtgGeneVector{x}=strArray{4};
                        
                    end
                    
                    if length(unique(qtgGeneVector))==1
                        
                        isQTGvector(m)=1;
                        geneVector{m}=qtgGeneVector{1};
                        
                        
                    else
                        
                        isQTNvector(m)=0;
                        isQTGvector(m)=0;
                        qtnPosVector(m)=0;
                        geneVector{m}=[];
                        varTypeVector{m}=[];
                        position{m}=[];
                        encoded{m}=[];  
                        chr{m}=[];
                        nucleotide(m)=0;
                        refCodon{m}=[];
                        altCodon{m}=[];
                        refDicodon1{m}=[];
                        refDicodon2{m}=[];
                        promoter{m}=[];
                        
                    end
                    
                else
                    
                    isQTNvector(m)=0;
                    isQTGvector(m)=0;
                    qtnPosVector(m)=0;
                    geneVector{m}=[];
                    varTypeVector{m}=[];
                    position{m}=[];
                    encoded{m}=[];  
                    chr{m}=[];
                    nucleotide(m)=0;
                    refCodon{m}=[];
                    altCodon{m}=[];
                    refDicodon1{m}=[];
                    refDicodon2{m}=[];
                    promoter{m}=[];
                    
                end
            else
                isQTNvector(m)=0;
                isQTGvector(m)=0;
                qtnPosVector(m)=0;
                geneVector{m}=[];
                varTypeVector{m}=[];
                position{m}=[];
                encoded{m}=[]; 
                chr{m}=[];
                nucleotide(m)=0; 
                refCodon{m}=[];
                altCodon{m}=[];
                refDicodon1{m}=[];
                refDicodon2{m}=[];
                promoter{m}=[];
            end
            
            
            
            m=m+1;
        end
        


    end
    
    
    
end





%now write tab-delimited output file
toOutput=table(conditionVector',timeVector',bPosVector',pValVector',betaVector',...
varExpVector',isQTNvector',qtnPosVector',isQTGvector',geneVector',varTypeVector',...
position',encoded',chr',nucleotide',refCodon',altCodon',...
promoter','VariableNames',...
{'condition','time','bPos','pVal','beta','varExp','isQTN','qtnPos','isQTG','gene','variantType',...
'genePosition','encoded','chr','nucleotide','refCodon','altCodon',...
'promoter'});
writetable(toOutput,'linearOutput.csv');








