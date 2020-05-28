function [OK_lambda,OK_tlag,OK_GamShape] = ABCrejection_tccDist_inference_BasalCellProlif(rtime,FreqDiv_experim,N,TolThres,lambda_avg,tlag_range,GamShape_range,xtickpos)
%% Approximate Bayesian Computation (ABC) algorithm to find best-fitting progenitor cell-cycle parameters for distributions in the No. of cell divisions
% Implementation of a single-rejection algorithm adapted to compare
% experimental H2BGFP-deconvoluted distributions of the No. of cell division
% rounds with inferred ones (coming from Monte-Carlo simulations of
% basal-cell turnover according to SP model).

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% FreqDiv_experim: {1,n} array containing experimental distributions in the No. of division rounds observed in basal cells at different time points (column vectors)
% N: No. of acceptable posterior estimates
% TolThres: Tolerance threshold for parameter value acceptance/rejection
% lambda_avg: average division rate (/week) (prior)
% tlag_range: range of possible values for the refractory period parameter (or minimum cell-cycle period, in weeks) to explore (prior)
% GamShape_range: range of possible values for the 'Shape' parameter of the gamma-distributed cell-cycle period to explore (prior)
% xtickpos: range of possible No. of division rounds

%% Output:
% OK_lambda: row vector of posterior, accepted values for the division rate parameter (/week)
% OK_tlag: row vector of posterior, accepted values for the refractory period length (weeks)
% OK_GamShape: row vector of posterior, accepted values for the 'Shape' parameter of the gamma-distributed cell-cycle period

%% Example:
% rtime = [0 1 2 3]; %(weeks)
% FreqDiv_experim{1,1} = hist(poissrnd(0,1000,1),[1:100])'; ...
%       FreqDiv_experim{1,2} = hist(poissrnd(1,1000,1),[1:100])'; ...
%       FreqDiv_experim{1,3} = hist(poissrnd(2,1000,1),[1:100])'; ...
%       FreqDiv_experim{1,4} = hist(poissrnd(3,1000,1),[1:100])';
% N = 100;
% TolThres = 0.15;
% lambda_avg = 1;
% tlag_range = [0:0.25:2]./7;
% GamShape_range = 2.^[0:6]; %[1,2,4,8,16,32,64];
% [OK_lambda,OK_tlag,OK_GamShape] = ABCrejection_tccDist_inference_BasalCellProlif(rtime,FreqDiv_experim,N,TolThres,lambda_avg,tlag_range,GamShape_range);

%% Initial definition of parameters:
M = 4000; % No. of simulated individual basal cells per trial
OK_lambda = zeros(1,N); OK_tlag = zeros(1,N); OK_GamShape = zeros(1,N);

%% ITERATION ON RANDOM PARAMETER VALUES TO FIND MATCHING SOLUTIONS
for aja = 1:N
    
    % Predefine distance metric:
    pdist = Inf;
    
    % Single-rejection method:
    while pdist > TolThres
        
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
        
    end
    
    % Record progression of the single-rejection iteration:
    disp(sprintf('%d/%d completed',aja,N))

    % Store the accepted parameter values (each estimate) in posterior distributions:
    OK_lambda(1,aja) = lambda;
    OK_tlag(1,aja) = tlag;
    OK_GamShape(1,aja) = GamShape;
     
end
