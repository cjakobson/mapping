%relies on output of parseVCF ('sgrpInputdata.mat')
clear all

load('sgrpInputData.mat')
clear condensedGenotype

clear coveredGenotype
clear coveredChroms

chromosomes=unique(inputMat.CHROM,'rows');

[nLoci,nStrains]=size(genotype);

strains={'BC187','DBVPG1106','DBVPG1373','DBVPG1788','DBVPG6044',...
    'DBVPG6765','L1374','L1528','SK1','UWOPS03-461.4',...
    'UWOPS83-787.3','UWOPS87-2421','W303','Y12','Y55',...
    'YJM975','YJM978','YJM981','YPS128'};


cleanStrains={'BC187','DBVPG6044','SK1','UWOPS03-461.4','UWOPS83-787.3','UWOPS87-2421',...
    'W303','Y12','Y55','YPS128','wine/European'};
%cleanStrains={'BC187','DBVPG6044','UWOPS03-461.4','Y12','YPS128','wine/European'};
%cleanStrains={'BC187','DBVPG6044','UWOPS03-461.4','UWOPS83-787.3','Y12','YPS128','wine/European'};
wineEuroIdx=[2 3 4 6 7 8 16 17 18]; %use these idx to condense wine/Euro genotypes (modal)
cleanIdx=[1 5 9 10 11 12 13 14 15 19];
%cleanIdx=[1 5 10 11 14 19];

cleanGenotype=genotype(:,cleanIdx);
wineEuroGenotype=mode(genotype(:,wineEuroIdx),2);

condensedGenotype=[cleanGenotype wineEuroGenotype];

[nLoci,nCondensedStrains]=size(condensedGenotype);

%remove loci with missing genotypes

coveredGenotype=condensedGenotype(sum(isnan(condensedGenotype),2)==0,:);
coveredChroms=inputMat.CHROM(sum(isnan(condensedGenotype),2)==0,:);
coveredPos=inputMat.POS(sum(isnan(condensedGenotype),2)==0,:);
coveredMutTypes=mutType(sum(isnan(condensedGenotype),2)==0,:);

[nLoci,nCoveredStrains]=size(coveredGenotype);

clear strainDist
clear nAlt
clear nInferred
m=1;
tic

window=250;
for i=(window+1):(nLoci-window)
    
    %idxToUse=ismember(coveredChroms(1:end,:),chromosomes(i,:),'rows');
    %build tree based on sliding window
    idxToUse=(i-window):(i+window);
    
    for j=1:nCoveredStrains
        for k=1:nCoveredStrains
        
            strainDist(j,k)=sum((coveredGenotype(idxToUse,j)-coveredGenotype(idxToUse,k)).^2,'omitnan');
        
        end
    end
    
    tree=seqneighjoin(strainDist,'equivar',cleanStrains);
    %subplot(4,4,i)
    if (i==7500)||(i==17500)
        %plot(tree,'orient','top')
        plot(tree,'Type','equalangle')
        title(['chromosome ' coveredChroms(i,:) 'position ' num2str(coveredPos(i))]);
    end
    
    
    genotypeToUse=coveredGenotype(idxToUse,:);
    
    %analyze inference of multiple parallel alleles using trees determined on a per-chromosome basis
    
    vInferred=inferParallel(genotypeToUse(window+1,:),cleanStrains,tree);
    nAlt(m)=sum(vInferred(1:length(cleanStrains))>0);
    nInferred(m)=sum(isnan(vInferred((length(cleanStrains)+1):end)));
    m=m+1;

    
end
toc

%generate chromosome II seq distance plot
idxToUse=ismember(coveredChroms(1:end,:),chromosomes(2,:),'rows');

%first compare Y55 to clean lineages
compareIdx=[2 4 8 10 11];
compareStrains=cleanStrains(compareIdx);

condensedGenotype=coveredGenotype(idxToUse,compareIdx);    
condensedPos=coveredPos(idxToUse);    


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
vY55=coveredGenotype(idxToUse,9);
hold on
for i=[1 5]
    plot(condensedPos,smooth(1-abs(condensedGenotype(:,i)-vY55),1000))
end
legend(compareStrains([1 5]))
ylim([0.5 1.2])
xlim([0 8*10^5])
title('Y55')

subplot(2,1,2)
vSK1=coveredGenotype(idxToUse,3);
hold on
for i=[1 3]
    plot(condensedPos,smooth(1-abs(condensedGenotype(:,i)-vSK1),1000))
end
legend(compareStrains([1 3]))
ylim([0.5 1.2])
xlim([0 8*10^5])
title('SK1')

clear nUnique
clear nShared

clear nParallel
clear nNotParallel

%now neutral expectation
for l=1:10
    tic
    load(['neutralAllGenotype' num2str(l) '.mat'])

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


figure('units','normalized','outerposition',[0 0 1 1])
%plot histograms
subplot(2,3,1)
histogram(nAlt)
title('number of strains with alt allele')
axis square

subplot(2,3,2)
histogram(nInferred(logical((nAlt>1).*(nAlt<nCondensedStrains))))
title('number of emergence events')
xlim([0 4])
axis square

%plot bars
subplot(2,3,4)
%bar([sum(nAlt==1)/length(nAlt),sum(nAlt~=1)/length(nAlt)])
bar([sum(nAlt==1),mean(nUnique);sum(nAlt~=1),mean(nShared)])
hold on
scatter(ones(1,10)*1.05+rand(1,10)*0.2,nUnique,'k','filled')
scatter(ones(1,10)*2.05+rand(1,10)*0.2,nShared,'k','filled')
xlim([0.5 2.5])
title('unique/shared')

subplot(2,3,5)
%bar([sum((nInferred==1).*(nAlt>1))/sum((nInferred>=1).*(nAlt>1)),sum((nInferred>1).*(nAlt>1))/sum((nInferred>=1).*(nAlt>1))])
bar([sum((nInferred==1).*(nAlt>1)),mean(nNotParallel);sum((nInferred>1).*(nAlt>1)),mean(nParallel)])
hold on
scatter(ones(1,10)*1.05+rand(1,10)*0.2,nNotParallel,'k','filled')
scatter(ones(1,10)*2.05+rand(1,10)*0.2,nParallel,'k','filled')
xlim([0.5 2.5])
title('emerging once/multiple')




%sliding analysis of variant type
vAll=variantTypes(coveredMutTypes);
vUnique=variantTypes(coveredMutTypes(nAlt==1));
vShared=variantTypes(coveredMutTypes(nAlt~=1));
vParallel=variantTypes(coveredMutTypes(nInferred>1));

subplot(2,3,6)
bar([vAll;vUnique;vShared;vParallel]')
xticklabels({'missense_variant','synonymous_variant','upstream_gene_variant',...
    'downstream_gene_variant','intergenic_region','other'})
 xtickangle(45)

 
 
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
plot(smooth(nAlt,10000))
hold on
plot(smooth(nAltNeutral,10000))
xlim([0 2.1*10^5])
ylim([1 4])



%output chromosome, pos, nAlt, and nEmerge
toOutput=table(coveredChroms((window+1):(end-window),:),...
    coveredPos((window+1):(end-window)),nAlt',nInferred','VariableNames',...
    {'chr','pos','nAlt','nInferred'});

writetable(toOutput,'sgrpInference.csv');




%find example of multiply emerging allele
[~, idx]=sort(nInferred);

%shift idx by window
idx=idx+250;

query=idx(end-100);

coveredChroms(query)
coveredPos(query)
coveredGenotype(query,:)
cleanStrains(coveredGenotype(query,:)==0)
cleanStrains(coveredGenotype(query,:)==1)


%plot local tree
idxToUse=(query-window):(query+window);
    
for j=1:nCoveredStrains
    for k=1:nCoveredStrains

        strainDist(j,k)=sum((coveredGenotype(idxToUse,j)-coveredGenotype(idxToUse,k)).^2,'omitnan');

    end
end

tree=seqneighjoin(strainDist,'equivar',cleanStrains);
plot(tree,'Type','equalangle')



figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
plot(smooth(nInferred,10000))
hold on
plot(smooth(nInferredNeutral,10000))
xlim([0 2.1*10^5])
ylim([1 1.5])

