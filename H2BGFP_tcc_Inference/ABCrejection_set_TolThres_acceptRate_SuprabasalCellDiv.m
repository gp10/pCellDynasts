function [pdist_all] = ABCrejection_set_TolThres_acceptRate_SuprabasalCellDiv(rtime,FreqDiv_experim,N,lambda_avg,tlagDiv,GamShapeDiv,mu_range,tlagShed_range,GamShapeShed_range,xtickpos)
%% Algorithm used to test distance metric values between simulated and experimental patterns in the number of cell division rounds in spinous (suprabasal) cells
% This is typically run as a pilot test before the Approximate Bayesian
% Computation (ABC) algorithm, so as to define a tolerance threshold for
% the ABC that sets an adequate level of acceptance rate.
% Simulations of distributions of no. of cell divisions in spinous layer are
% computed similarly as in the ABC algorithm (coming from Monte-Carlo
% simulations of suprabasal cell replacement by basal cycling progenitors
% according to SP model) and the same distance metric is preserved.

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% FreqDiv_experim: {1,n} array containing experimental distributions in the No. of division rounds observed in spinous cells at different time points (column vectors)
% N: No. of trials
% lambda_avg: average division rate (/week) (prior)
% tlagDiv: refractory period (or minimum cell-cycle period, in weeks) before division (prior)
% GamShapeDiv: 'Shape' parameter of the gamma-distributed cell-cycle period for division (prior)
% mu_range: range of possible values for average shedding/terminal differentiation rate of spinous cells (/week) (prior)
% tlagShed_range: range of possible refractory periods (or minimum spinous-layer residency periods) before shedding/terminal differentiation, expressed relative to avg. stratification/shedding time (1/mu) (prior)
% GamShapeShed_range: range of possible values for the 'Shape' parameter of the gamma-distributed residency periods in spinous layer (prior)
% xtickpos: range of possible No. of division rounds

%% Output:
% pdist_all: row vector of distance metric values of experimental vs. simulated frequencies in the No. of cell division rounds in spinous compartment under random waiting-time parameter values

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
M = 400; % No. of simulated individual suprabasal cells per trial
pdist_all = zeros(1,N);

%% ITERATION ON RANDOM PARAMETER VALUES TO RETRIEVE GOODNESS-OF-FIT DISTANCE METRIC
for aja = 1:N
    
    % Predefine distance metric:
    pdist = Inf;
    
    % Draw random parameter values from prior (discrete) distributions:
    lambda = lambda_avg; %2.8+rand*1.2;
    tlagDiv;
    GamShapeDiv;
    mu = mu_range(find(mnrnd(1,1/length(mu_range)*ones(1,length(mu_range)))==1));
    tlagShed = (rand*(tlagShed_range(end)-tlagShed_range(1)) + tlagShed_range(1)).*1/mu; % e.g. up to 2d maximum
    GamShapeShed = GamShapeShed_range(find(mnrnd(1,1/length(GamShapeShed_range)*ones(1,length(GamShapeShed_range)))==1));
    rho = mu / (lambda + mu); % homeostatic requirement

    % Simulation of SP-model time course in the No. of division rounds in spinous (suprabasal) compartment:
    [NoDiv_SL] = MonteCarloSimulator_SP_SuprabasalCellDiv(rtime,M,rho,lambda,tlagDiv,GamShapeDiv,mu,tlagShed,GamShapeShed);

    % actual division values are translated into freqs (FreqDiv_SPmodel) to be compared against experimental freqs (FreqDiv_experim)
    if iscell(xtickpos) == 0 % if the max. No of divisions tracked is the same for all time points
        myspan = [xtickpos-1 xtickpos(end)];
        for buc = 1:length(rtime)
            [counts,centres] = hist(NoDiv_SL(buc,:),myspan);
            FreqDiv_SPmodel{1,buc} = counts(1:end-1) ./ sum(counts(1:end-1));
        end
    else % max. No of divisions tracked is different for each time point
        for buc = 1:length(rtime)
            myspan = [xtickpos{1,buc}-1 xtickpos{1,buc}(end)];
            [counts,centres] = hist(NoDiv_SL(buc,:),myspan);
            FreqDiv_SPmodel{1,buc} = counts(1:end-1) ./ sum(counts(1:end-1));
        end
    end

    % Distance metric of simulated vs. experimental histograms in the No. of cell division rounds in spinous layer (based on Residual Sum of Squares (RSS)):
    pdist_time = [];
    for buc = 1:length(rtime)
        pdist_time(buc,1) = sum((FreqDiv_SPmodel{1,buc} - FreqDiv_experim{1,buc}').^2);
    end
    pdist = sum(pdist_time,1);
    
    % Store distance metric value
    pdist_all(1,aja) = pdist;
    
end
