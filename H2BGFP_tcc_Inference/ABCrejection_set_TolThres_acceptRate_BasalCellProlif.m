function [pdist_all] = ABCrejection_set_TolThres_acceptRate_BasalCellProlif(rtime,FreqDiv_experim,N,lambda_avg,tlag_range,GamShape_range,xtickpos)
%% Algorithm used to test distance metric values between simulated and experimental patterns in the number of cell division rounds
% This is typically run as a pilot test before the Approximate Bayesian
% Computation (ABC) algorithm, so as to define a tolerance threshold for
% the ABC that sets an adequate level of acceptance rate.
% Simulations of distributions of no. of cell divisions are computed
% similarly as in the ABC algorithm (coming from Monte-Carlo simulations of
% basal-cell turnover according to SP model) and the same distance metric is
% preserved.

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% FreqDiv_experim: {1,n} array containing experimental distributions in the No. of division rounds observed in basal cells at different time points (column vectors)
% N: No. of trials (random tcc parameter combinations)
% lambda_avg: average division rate (/week) (prior)
% tlag_range: range of possible values for the refractory period parameter (or minimum cell-cycle period, in weeks) to explore (prior)
% GamShape_range: range of possible values for the 'Shape' parameter of the gamma-distributed cell-cycle period to explore (prior)
% xtickpos: range of possible No. of division rounds

%% Output:
% pdist_all: row vector of distance metric values of experimental vs. simulated frequencies in the No. of cell division rounds under random tcc parameter values

%% Example:
% rtime = [0 1 2 3]; %(weeks)
% FreqDiv_experim{1,1} = hist(poissrnd(0,1000,1),[1:100])'; ...
%       FreqDiv_experim{1,2} = hist(poissrnd(1,1000,1),[1:100])'; ...
%       FreqDiv_experim{1,3} = hist(poissrnd(2,1000,1),[1:100])'; ...
%       FreqDiv_experim{1,4} = hist(poissrnd(3,1000,1),[1:100])';
% N = 100;
% lambda_avg = 1;
% tlag_range = [0:0.25:2]./7;
% GamShape_range = 2.^[0:6]; %[1,2,4,8,16,32,64];
% [pdist_all] = ABCrejection_set_TolThres_acceptRate_BasalCellProlif(rtime,FreqDiv_experim,N,lambda_avg,tlag_range,GamShape_range)

%% Initial definition of parameters:
M = 4000; % No. of simulated individual basal cells per trial
pdist_all = zeros(1,N);

%% ITERATION ON RANDOM PARAMETER VALUES TO RETRIEVE GOODNESS-OF-FIT DISTANCE METRIC
for aja = 1:N
    
    % Predefine distance metric:
    pdist = Inf;
    
    % Draw random parameter values from prior (discrete) distributions:
    lambda = lambda_avg; %2.8+rand*1.2;
    tlag = rand*(tlag_range(end)-tlag_range(1)) + tlag_range(1); % e.g. up to 2d maximum
    GamShape = GamShape_range(find(mnrnd(1,1/length(GamShape_range)*ones(1,length(GamShape_range)))==1));

    % Simulation of SP-model time course in the No. of division rounds:
    [NoDiv_BL] = MonteCarloSimulator_SP_BasalCellProlif(rtime,lambda,M,tlag,GamShape);

    % actual division values are translated into freqs (FreqDiv_SPmodel) to be compared against experimental freqs (FreqDiv_experim)
    if iscell(xtickpos) == 0 % if the max. No of divisions tracked is the same for all time points
        myspan = [xtickpos-1 xtickpos(end)];
        for buc = 1:length(rtime)
            [counts,centres] = hist(NoDiv_BL(buc,:),myspan);
            FreqDiv_SPmodel{1,buc} = counts(1:end-1) ./ sum(counts(1:end-1));
        end
    else % max. No of divisions tracked is different for each time point
        for buc = 1:length(rtime)
            myspan = [xtickpos{1,buc}-1 xtickpos{1,buc}(end)];
            [counts,centres] = hist(NoDiv_BL(buc,:),myspan);
            FreqDiv_SPmodel{1,buc} = counts(1:end-1) ./ sum(counts(1:end-1));
        end
    end

    % Distance metric of simulated vs. experimental histograms in the No. of cell division rounds (based on Residual Sum of Squares (RSS)):
    pdist_time = [];
    for buc = 1:length(rtime)
        pdist_time(buc,1) = sum((FreqDiv_SPmodel{1,buc} - FreqDiv_experim{1,buc}').^2);
    end
    pdist = sum(pdist_time,1);
    
    % Store distance metric value
    pdist_all(1,aja) = pdist;
    
end
