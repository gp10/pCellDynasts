function [OK_mu,OK_tlagShed,OK_GamShapeShed] = ABCrejection_tccDist_inference_SuprabasalCellDiv(rtime,FreqDiv_experim,N,TolThres,lambda_avg,tlagDiv,GamShapeDiv,mu_range,tlagShed_range,GamShapeShed_range,xtickpos)
%% Approximate Bayesian Computation (ABC) algorithm to find best-fitting progenitor stratification/shedding waiting-time parameters for distributions in the No. of cell divisions
% Implementation of a single-rejection algorithm adapted to compare
% experimental H2BGFP-deconvoluted distributions of the No. of cell division
% rounds experienced by suprabasal cells with inferred ones (coming from
% Monte-Carlo simulations of suprabasal cell replacement by basal cycling
% progenitors according to SP model).

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% FreqDiv_experim: {1,n} array containing experimental distributions in the No. of division rounds observed in spinous cells at different time points (column vectors)
% N: No. of acceptable posterior estimates
% TolThres: Tolerance threshold for parameter value acceptance/rejection
% lambda_avg: average division rate (/week) (prior)
% tlagDiv: refractory period (or minimum cell-cycle period, in weeks) before division (prior)
% GamShapeDiv: 'Shape' parameter of the gamma-distributed cell-cycle period for division (prior)
% mu_range: range of possible values for average shedding/terminal differentiation rate of spinous cells (/week) (prior)
% tlagShed_range: range of possible refractory periods (or minimum spinous-layer residency periods) before shedding/terminal differentiation, expressed relative to avg. stratification/shedding time (1/mu) (prior)
% GamShapeShed_range: range of possible values for the 'Shape' parameter of the gamma-distributed residency periods in spinous layer (prior)
% xtickpos: range of possible No. of division rounds

%% Output:
% OK_mu: row vector of posterior, accepted values for the shedding rate parameter (/week)
% OK_tlagShed: row vector of posterior, accepted values for the refractory period length (weeks)
% OK_GamShapeShed: row vector of posterior, accepted values for the 'Shape' parameter of the gamma-distributed residency periods in spinous layer

%% Example:
% rtime = [0 1 2 3]; %(weeks)
% FreqDiv_experim{1,1} = hist(poissrnd(0,1000,1),[1:100])'; ...
%       FreqDiv_experim{1,2} = hist(poissrnd(1,1000,1),[1:100])'; ...
%       FreqDiv_experim{1,3} = hist(poissrnd(2,1000,1),[1:100])'; ...
%       FreqDiv_experim{1,4} = hist(poissrnd(3,1000,1),[1:100])';
% N = 100;
% TolThres = 0.15;
% lambda_avg = 1;
% tlagDiv = 0.5/7;
% GamShapeDiv = 1;
% mu_range = [2.^[0:1:5]].*lambda_avg;
% tlagShed_range = [0:0.3:0.6]; %given in relative terms to avg. stratification/shedding time (1/mu)
% GamShapeShed_range = 2.^[0:0.5:2];
% [OK_mu,OK_tlagShed,OK_GamShapeShed] = ABCrejection_tccDist_inference_SuprabasalCellDiv(rtime,FreqDiv_experim,N,TolThres,lambda_avg,tlagDiv,GamShapeDiv,mu_range,tlagShed_range,GamShapeShed_range,xtickpos);

%% Initial definition of parameters:
M = 400; % No. of simulated individual suprabasal cells per trial
OK_mu = zeros(1,N); OK_tlagShed = zeros(1,N); OK_GamShapeShed = zeros(1,N);

%% ITERATION ON RANDOM PARAMETER VALUES TO FIND MATCHING SOLUTIONS
for aja = 1:N
    
    % Predefine distance metric:
    pdist = Inf;
    
    % Single-rejection method:
    while pdist > TolThres
        
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
        
    end
    
    % Record progression of the single-rejection iteration:
    disp(sprintf('%d/%d completed',aja,N))

    % Store the accepted parameter values (each estimate) in posterior distributions:
    OK_mu(1,aja) = mu;
    OK_tlagShed(1,aja) = tlagShed;
    OK_GamShapeShed(1,aja) = GamShapeShed;
     
end
