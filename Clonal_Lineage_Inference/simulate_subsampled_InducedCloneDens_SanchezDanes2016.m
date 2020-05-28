function [sampled_avgNClones, sampled_stdNClones] = simulate_subsampled_InducedCloneDens_SanchezDanes2016(CloneSizes,NoMicePersis,NoClones_mean_t0,NoClones_std_t0)
%% RETRIEVES DESCRIPTIVE STATISTICS ON CLONE DENSITY FROM SIMULATIONS STARTING WITH A VARIABLE LABELLING INDUCTION:
% The experimental variability in the number of clones observed per
% individual mouse at the earliest time point is reproduced, and then
% time-course clone simulations dissected to retrieve the number of
% persisting clones for each simulated mouse at a given time
% post-induction, ultimately retrieving descriptive statistics on clonal
% density values across simulated mice and time.

% from Piedrafita et al, 2020

%% Input:
% CloneSizes: mxn matrix of clone sizes, containing the number of basal cells of m clones at n time points
% NoMicePersis: row vector containing the number of mice from which to extract statistics at each time point (same as experimental number of mice)
% NoClones_mean_t0: column vector containing the average number of experimental clones counted per unit area per animal at the first time point (t0)
% NoClones_std_t0: column vector containing the standard deviation of the number of experimental clones counted per unit area per animal at the first time point (t0)

%% Output:
% sampled_avgNClones: row vector containing the mean number of persisting clones across animals at the different time points
% sampled_stdNClones: row vector containing the standard deviation of the number of persisting clones across animals at the different time points

%%
sampled_NClones = {};
sampled_avgNClones = []; sampled_stdNClones = [];
ScalingFactor_t0 = size(CloneSizes,1) / size(find(CloneSizes(:,1)~=0),1); % to correct for clones lost by first experimental time point (t0)
% Independent random sampling for each mouse at each time point:
for aja = 1:size(CloneSizes,2)
    for eje = 1:NoMicePersis(aja)
        % simulate a random initial label induction efficiency:
        label_induc_t0 = round((NoClones_std_t0*randn+NoClones_mean_t0)*ScalingFactor_t0);
        % simulate the time course in the number of persisting clones starting from a random set of as many clones as indicated by the induction
        mypermut = randperm(size(CloneSizes,1));
        sampled_NClones{1,aja}(eje,1) = size(find(CloneSizes(mypermut(1:label_induc_t0),aja)~=0),1);
    end
    % Retrieve descriptive statistics (avg number of clones and std):
    sampled_avgNClones(1,aja) = mean(sampled_NClones{1,aja});
    sampled_stdNClones(1,aja) = std(sampled_NClones{1,aja},0);
end
