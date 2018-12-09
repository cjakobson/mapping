%relies on output of parseVCF ('sgrpInputdata.mat')

clear all


load('sgrpInputData.mat')
clear condensedGenotype

clear coveredParadoxusGenotype
clear coveredChroms

[nLoci,nStrains]=size(paradoxusGenotype);

paradoxusStrains={'IFO1804', 'KPN3829', 'N0x2D44', 'Q310x2E4',...
    'Q320x2E3', 'Q590x2E1', 'Q620x2E5', 'Q690x2E8', 'Q740x2E4',...
    'Q890x2E8', 'Q950x2E3', 'S360x2E7', 'T210x2E4', 'UFRJ50816',...
    'W7', 'Y60x2E5', 'Y7', 'Y80x2E1', 'Y8_5', 'Y9_6', 'YPS138',...
    'Z1', 'Z1_1'};


%remove loci with missing genotypes
coveredParadoxusGenotype=paradoxusGenotype(sum(isnan(paradoxusGenotype),2)<3,:);
coveredParadoxusChroms=inputMat3.CHROM(sum(isnan(paradoxusGenotype),2)<3,:);

[nLoci,nCoveredStrains]=size(coveredParadoxusGenotype);

window=250;

clear strainDist
clear nAlt
clear nInferred
m=1;
tic
for i=(window+1):(nLoci-window)
    
    %idxToUse=ismember(coveredChroms(1:end,:),chromosomes(i,:),'rows');
    %build tree based on 500 variant sliding window
    idxToUse=(i-window):(i+window);
    
    for j=1:nCoveredStrains
        for k=1:nCoveredStrains
        
            strainDist(j,k)=sum((coveredParadoxusGenotype(idxToUse,j)-coveredParadoxusGenotype(idxToUse,k)).^2,'omitnan');
        
        end
    end
    
    tree=seqneighjoin(strainDist,'equivar',paradoxusStrains);
    %subplot(4,4,i)
    %plot(tree,'orient','top')
    %plot(tree,'Type','equalangle')
    %title(['chromosome ' chromosomes(i,:)]);
    
    genotypeToUse=coveredParadoxusGenotype(idxToUse,:);
    
    %analyze inference of multiple parallel alleles using trees determined on a per-chromosome basis
    
    vInferred=inferParallel(genotypeToUse(101,:),paradoxusStrains,tree);
    nAlt(m)=sum(vInferred(1:length(paradoxusStrains))>0);
    nInferred(m)=sum(isnan(vInferred((length(paradoxusStrains)+1):end)));
    m=m+1;
    if mod(m,10000)==0
        m
    end

    
end
toc


clear nUnique
clear nShared

clear nParallel
clear nNotParallel

%now neutral expectation
for l=1:10
    tic
    load(['neutralParadoxusGenotype' num2str(l) '.mat'])

    [nNeutralLoci,nNeutralStrains]=size(neutralGenotype);
    
    nAltNeutral=zeros(nNeutralLoci,1);
    nInferredNeutral=zeros(nNeutralLoci,1);
    
    for j=1:nNeutralStrains
        for k=1:nNeutralStrains
        
            strainDist(j,k)=sum((neutralGenotype(:,j)-neutralGenotype(:,k)).^2,'omitnan');
        
        end
    end
    
    neutralTree=seqneighjoin(strainDist,'equivar');
    neutralNames=get(neutralTree,'NodeNames');
    
    clear strainDist
    clear vInferred
    m=1;
    
    for i=1:nNeutralLoci
        %analyze inference of multiple parallel alleles 

        vInferred=inferParallel(neutralGenotype(i,:),neutralNames,neutralTree);
        nAltNeutral(m)=sum(vInferred(1:nNeutralStrains)>0);
        nInferredNeutral(m)=sum(isnan(vInferred((nNeutralStrains+1):end)));
        m=m+1;

    
    end
    


    nUnique(l)=sum(nAltNeutral==1);
    nShared(l)=sum(nAltNeutral~=1);
    
    nParallel(l)=sum((nInferredNeutral>1).*(nAltNeutral>1));
    nNotParallel(l)=sum((nInferredNeutral==1).*(nAltNeutral>1));
    
    toc
end



figure
%plot histograms
subplot(2,3,1)
histogram(nAlt)
title('number of strains with alt allele')

subplot(2,3,2)
histogram(nInferred(logical((nAlt>1).*(nAlt<nCoveredStrains))))
title('number of emergence events')
xlim([0 4])

%plot bars
subplot(2,3,4)
%bar([sum(nAlt==1)/length(nAlt),sum(nAlt~=1)/length(nAlt)])
bar([sum(nAlt==1),mean(nUnique);sum(nAlt~=1),mean(nShared)])
hold on
scatter(ones(1,10)*1.05+rand(1,10)*0.2,nUnique,'k','filled')
scatter(ones(1,10)*2.05+rand(1,10)*0.2,nShared,'k','filled')
%ylim([0 1])
xlim([0.5 2.5])
title('unique/shared')

subplot(2,3,5)
%bar([sum((nInferred==1).*(nAlt>1))/sum((nInferred>=1).*(nAlt>1)),sum((nInferred>1).*(nAlt>1))/sum((nInferred>=1).*(nAlt>1))])
bar([sum((nInferred==1).*(nAlt>1)),mean(nNotParallel);sum((nInferred>1).*(nAlt>1)),mean(nParallel)])
hold on
scatter(ones(1,10)*1.05+rand(1,10)*0.2,nNotParallel,'k','filled')
scatter(ones(1,10)*2.05+rand(1,10)*0.2,nParallel,'k','filled')
%ylim([0 1])
xlim([0.5 2.5])
title('emerging once/multiple')

