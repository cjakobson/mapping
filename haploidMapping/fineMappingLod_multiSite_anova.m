function [p_anova,pairwise_p,ph1,ph2] = fineMappingLod_multiSite_anova(pos1,pos2,lower,upper,queryPos,x,y,b_fwselection,residual)
yr = residual;
yr = yr + b_fwselection(queryPos)*x(:,queryPos);

groupNames = {'YJM/YJM','YJM/RM','RM/YJM','RM/RM'};
groupAssignments{1} = x(:,pos1)==-1 & x(:,pos2) == -1;
groupAssignments{2} = x(:,pos1)==-1 & x(:,pos2) == 1;
groupAssignments{3} = x(:,pos1)==1 & x(:,pos2) == -1;
groupAssignments{4} = x(:,pos1)==1 & x(:,pos2) == 1;

group = cell(length(yr),1);
for i = 1:length(yr)
    for j = 1:4
        if groupAssignments{j}(i) == 1
            group{i} = groupNames{j};
        end
    end
end
toRemove = find(cellfun(@isempty,group));
newx = x;
newx(toRemove,:) = [];
yr(toRemove) = [];
group(toRemove) = [];
for i = 1:4
    groupAssignments{i}(toRemove) = [];
end

%%% Do not run the ANOVA analysis for low coverage genotypes with < 200 of
%%%  RM/RM and YJM/YJM genotypes, or those with no crossover genotypes.
%%%  Requires both types of crossovers.
if sum(groupAssignments{1}) > length(y)/4 && sum(groupAssignments{4}) > length(y)/4
    if sum(groupAssignments{2}) > 0 && sum(groupAssignments{3}) > 0
        
        % Use anova on the 4 genotype groups at position 1 and position 2
        [p_anova,tbl_anova,stats_anova] = anova1(yr,group,'off');
        
        %%% Perform Anova comparison of equality for all pairs of genotypes.
        [pairwise_p,mean_group,handle,gnames] = multcompare(stats_anova,'Display','off');
        
        %%%% Group order in pairwise_p is messed up in unpredictable ways -
        %%%% but at least it's outputted to gnames. Make an indexing
        %%%% variable.
        index1 = find(strcmp(groupNames{1},gnames));
        index2 = find(strcmp(groupNames{2},gnames));
        index3 = find(strcmp(groupNames{3},gnames));
        index4 = find(strcmp(groupNames{4},gnames));
        %%% Index order for first two columns of pairwise_p is lower/upper
        %%% always. But we don't know if index1 > index2 or vice versa.
        ph1_index1 = find(pairwise_p(:,1) == min(index1,index2) & pairwise_p(:,2) == max(index1,index2));
        ph1_index2 = find(pairwise_p(:,1) == min(index3,index4) & pairwise_p(:,2) == max(index3,index4));
        ph2_index1 = find(pairwise_p(:,1) == min(index1,index3) & pairwise_p(:,2) == max(index3,index1));
        ph2_index2 = find(pairwise_p(:,1) == min(index2,index4) & pairwise_p(:,2) == max(index4,index2));
        
        %%% ph1 stands for p-value H1, the likelihood that hypothesis 1 is
        %%% true. H1 is that pos1 is the causal variant (i.e. YJM/YJM ==
        %%% YJM/RM and RM/RM == RM/YJM)
        ph1 = pairwise_p(ph1_index1,6)*pairwise_p(ph1_index2,6);
        
        %%% ph2 stands for p-value H2, the likelihood that hypothesis 2 is
        %%% true. H1 is that pos2 is the causal variant. This is the alternate hypothesis. (i.e. YJM/YJM ==
        %%% RM/YJM and RM/RM == YJM/RM)
        ph2 = pairwise_p(ph2_index1,6)*pairwise_p(ph2_index2,6);
    else
        p_anova = -1;
        pairwise_p = -1;
        ph1 = -1;
        ph2 = -1;
    end
else
    p_anova = -1;
    pairwise_p = -1;
    ph1 = -1;
    ph2 = -1;
end
