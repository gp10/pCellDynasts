function [x,tini,tend,CloneSizePicked] = sampling_Initial_Clone_Fullgrabe2015(iniCloneSizes,dens)
%% RANDOM SAMPLING OF AN INITIAL CLONE SIZE AND COMPOSITION TO BE SIMULATED:
% A preformed clone is picked with a number of labelled cells consistent
% with clone size frequencies at 40d post-induction in Fullgrabe et al's
% experiments. A random cellular composition is assigned consistent with
% the global proportion of progenitor cells in the basal layer (rho).

% from Piedrafita et al, 2020

%% Input:
% iniCloneSizes: vector of [NClones,1] shape containing initial experimental clone sizes
% dens: numeric; global proportion of progenitor cells in the basal layer (rho)

%% Output:
% x: matrix of [CloneSizePicked,2] shape specifying P,D type of individual cells within the clone
% tini: default value for the time of birth of constituent cells
% tend: default value for the time of death of constituent cells
% CloneSizePicked: number of basal cells assigned to the preformed clone

%%
% Clone Freqs at rtime(3)=5.7143w (used to estimate initial clone size) are retreived:
myfreq = histc(iniCloneSizes,[1:1:max(max(iniCloneSizes))],1);
mfreq_rel = myfreq ./ sum(myfreq,1);
sem_mfreq_rel = sqrt(mfreq_rel .* (1-mfreq_rel)) ./ sqrt(sum(myfreq,1));

% Pick a random clone size:
CloneSize_Prob = [];
for eje = 1:size(mfreq_rel,1)
    CloneSize_Prob(eje) = sem_mfreq_rel(eje,1).*randn + mfreq_rel(eje,1);
    if CloneSize_Prob(eje)<0
        CloneSize_Prob(eje) = 0;
    end
end
CloneSize_Prob_norm = CloneSize_Prob ./ sum(CloneSize_Prob);
CloneSizePicked = find(mnrnd(1,CloneSize_Prob_norm));

% Pick a random clone composition:
x = [];
for elem = 1:CloneSizePicked
    if dens >= rand
        x = [x; 1 0];
    else
        x = [x; 0 1];
    end
end

% Predefine birth and death time variables for individual cells
tini = zeros(CloneSizePicked,1);
tend = zeros(CloneSizePicked,1);
