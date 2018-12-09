%function that returns variant type distribution
function [distOut,arrayOut]=variantTypes(inputArray)

variantArray={'missense_variant','synonymous_variant','upstream_gene_variant',...
    'downstream_gene_variant','intergenic_region'};

arrayOut=zeros(1,length(inputArray));

for i=1:length(inputArray)

    toParse=inputArray(i);
    
    for j=1:length(variantArray)
        
        if strcmp(toParse,variantArray(j))
            arrayOut(i)=j;
        end
        
    end

end

arrayOut(arrayOut==0)=length(variantArray)+1;

for i=1:length(variantArray)+1
    
    distOut(i)=sum(arrayOut==i);
    
end

distOut=distOut./sum(distOut);

