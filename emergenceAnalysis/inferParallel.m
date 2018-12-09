%function to infer parallel allele occurrence from phylo tree and
%genotypes

function vInferred=inferParallel(vInputGenotype,vInputStrains,inputTree)

%first retrieve some info about the phylogeny

children=get(inputTree,'Pointers'); %happily these are listed leaf-to-root, so no need to reorder
nodes=get(inputTree,'NodeNames');   %this includes leaves and branches
branches=get(inputTree,'BranchNames');
leaves=get(inputTree,'LeafNames');

for i=1:length(leaves)  %reindex to convert input genotypes into order used in tree
    leafIdx(i)=find(ismember(vInputStrains,leaves(i)));
    vInferred(i)=vInputGenotype(leafIdx(i));
end

for i=1:length(branches)    %impute genotypes/conficts at branch points
    
    childIdx=children(i,:);
    
    if vInferred(childIdx(1))==vInferred(childIdx(2))   %keep consensus genotype
        vInferred(length(leaves)+i)=vInferred(childIdx(1));
    elseif isnan(vInferred(childIdx(1)))&&isnan(vInferred(childIdx(2))) %if both are clash, register mode
        vInferred(length(leaves)+i)=mode(vInputGenotype);
    elseif isnan(vInferred(childIdx(1))) %if only one is a clash, keep the other
        vInferred(length(leaves)+i)=vInferred(childIdx(2));
    elseif isnan(vInferred(childIdx(2))) %if only the other is a clash, keep the first
        vInferred(length(leaves)+i)=vInferred(childIdx(1));
    elseif vInferred(childIdx(1))~=vInferred(childIdx(2))   %if unequal, register clash
        vInferred(length(leaves)+i)=NaN;
    else    %otherwise revert to ancestral
        vInferred(length(leaves)+i)=0;
    end

end


end