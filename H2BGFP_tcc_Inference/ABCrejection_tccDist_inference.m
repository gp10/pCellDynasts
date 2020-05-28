function [OK_lambda,OK_tlag,OK_GamShape] = ABCrejection_tccDist_inference(rtime,rnlog2_H2BGFP_all,N,TolThres,lambda_avg,tlag_range,GamShape_range)
%% Approximate Bayesian Computation (ABC) algorithm to find best-fitting progenitor cell-cycle parameters for H2BGFP dilution pattern
% Implementation of a single-rejection algorithm adapted to compare
% experimental H2BGFP intensity distributions with inferred ones (coming
% from Monte-Carlo simulations of basal-cell turnover according to SP model).

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% rnlog2_H2BGFP_all: {1,n} array containing experimental distributions of H2BGFP intensities for the different time points (expressed as log2 values; column vectors)
% N: No. of acceptable posterior estimates
% TolThres: Tolerance threshold for parameter value acceptance/rejection
% lambda_avg: average division rate (/week) (prior)
% tlag_range: range of possible values for the refractory period parameter (or minimum cell-cycle period, in weeks) to explore (prior)
% GamShape_range: range of possible values for the 'Shape' parameter of the gamma-distributed cell-cycle period to explore (prior)

%% Output:
% OK_lambda: row vector of posterior, accepted values for the division rate parameter (/week)
% OK_tlag: row vector of posterior, accepted values for the refractory period length (weeks)
% OK_GamShape: row vector of posterior, accepted values for the 'Shape' parameter of the gamma-distributed cell-cycle period

%% Example:
% rtime = [0 1 2 3]; %(weeks)
% rnlog2_H2BGFP_all{1,1} = normrnd(0,0.2,1000,1); ...
%       rnlog2_H2BGFP_all{1,2} = normrnd(-1,0.2,1000,1); ...
%       rnlog2_H2BGFP_all{1,3} = normrnd(-2,0.2,1000,1); ...
%       rnlog2_H2BGFP_all{1,4} = normrnd(-3,0.2,1000,1);
% N = 100;
% TolThres = 0.15;
% lambda_avg = 1;
% tlag_range = [0:0.25:2]./7;
% GamShape_range = 2.^[0:6]; %[1,2,4,8,16,32,64];
% [OK_lambda,OK_tlag,OK_GamShape] = ABCrejection_tccDist(rtime,rnlog2_H2BGFP_all,N,TolThres,lambda_avg,tlag_range,GamShape_range);

%% Initial definition of parameters:
M = 1000; % No. of simulated individual basal cells per trial
OK_lambda = zeros(1,N); OK_tlag = zeros(1,N); OK_GamShape = zeros(1,N);

% Experimental H2BGFP intensity values are normalized to median at time0
rlog2_Idist = {};
for buc = 1:length(rtime)
    rlog2_Idist{1,buc} = rnlog2_H2BGFP_all{1,buc}-median(rnlog2_H2BGFP_all{1,1},1);
end

%% ITERATION ON RANDOM PARAMETER VALUES TO FIND MATCHING SOLUTIONS
for aja = 1:N
    
    % Predefine distance metric:
    pdist = 0;
    
    % Single-rejection method:
    while pdist < TolThres
        
        % Draw random parameter values from prior (discrete) distributions:
        lambda = lambda_avg; %2.8+rand*1.2;
        tlag = rand*(tlag_range(end)-tlag_range(1)) + tlag_range(1); % e.g. up to 2d maximum
        GamShape = GamShape_range(find(mnrnd(1,1/length(GamShape_range)*ones(1,length(GamShape_range)))==1));

        % Simulation of H2BGFP dilution time course:
        Idist = MonteCarloSimulator_SP_BasalCellProlif_H2BGFPdil(rtime,rnlog2_H2BGFP_all{1,1},lambda,M,tlag,GamShape);

        % Simulated H2BGFP intensity values are translated into log2 values and these normalized to avg at time0 (as in experimental data)
        clog2_Idist = log2(Idist)-median(log2(Idist(1,:)),2);

        % Distance metric of simulated vs. experimental H2BGFP intensity distributions (based on Qs):
        pdist_time = [];
        for buc = 1:length(rtime)
            pdist_time(buc,1) = sum(abs(quantile(clog2_Idist(buc,:),[0.025 0.25 0.5 0.75 0.975])-quantile(rlog2_Idist{1,buc},[0.025 0.25 0.5 0.75 0.975])));
        end
        pdist = 1/sum(pdist_time,1);
        
    end
    
    % Record progression of the single-rejection iteration:
    disp(sprintf('%d/%d completed',aja,N))

    % Store the accepted parameter values (each estimate) in posterior distributions:
    OK_lambda(1,aja) = lambda;
    OK_tlag(1,aja) = tlag;
    OK_GamShape(1,aja) = GamShape;
     
end
