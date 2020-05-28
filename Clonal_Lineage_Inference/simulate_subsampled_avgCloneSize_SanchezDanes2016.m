function [sampled_avgCloneSize, sampled_semCloneSize] = simulate_subsampled_avgCloneSize_SanchezDanes2016(CloneSizes,NoTrackedClones)
%% RETRIEVES DESCRIPTIVE CLONE SIZE STATISTICS FROM A SPECIFIED NUMBER OF CLONES OBTAINED BY RANDOM SUBSAMPLING OF A LARGE POOL OF CLONAL DATA SETS:
% An average and sem values are calculated for random subsets of clones
% at any given time point, the number of subsampled clones being determined
% previously by the user.

% from Piedrafita et al, 2020

%% Input:
% CloneSizes: mxn matrix of clone sizes, containing the number of basal cells of m clones at n time points
% NoTrackedClones: row vector containing the number of random clones from which to extract statistics at each time point

%% Output:
% sampled_avgCloneSize: row vector containing the mean size from clones sampled at different time points
% sampled_semCloneSize: row vector containing the sem size from clones sampled at different time points

%%
% Independent random sampling for each time point:
sampled_CloneSizes = {};
for aja = 1:size(CloneSizes,2)
    mypermut = randperm(size(CloneSizes,1));
    rnd_pickpos = find(CloneSizes(mypermut,aja)~=0,NoTrackedClones(aja));
    sampled_CloneSizes{1,aja} = CloneSizes(mypermut(rnd_pickpos),aja);
end
% Retrieve descriptive statistics (avg sample clone size and sem):
for ete = 1:size(sampled_CloneSizes,2)
    sampled_avgCloneSize(ete) = mean(sampled_CloneSizes{1,ete});
    sampled_semCloneSize(ete) = std(sampled_CloneSizes{1,ete},0) ./ sqrt(size(sampled_CloneSizes{1,ete},1));
end
