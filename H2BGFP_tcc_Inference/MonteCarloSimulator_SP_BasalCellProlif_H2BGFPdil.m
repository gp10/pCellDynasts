function [Idist] = MonteCarloSimulator_SP_BasalCellProlif_H2BGFPdil(rtime,It0,lambda,indiv,tlag,GamShape)
%% Non-Markovian Monte Carlo simulator of H2BGFP intensity dilution in basal cells
% H2BGFP intensity is tracked over time as an initially labelled
% progenitor sucessively divides (only a one-branch pathway -the intensity
% of a given sibling- is considered upon each division).
% It allows to consider gamma-distributed cell cycle periods.

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks) - the 1st is taken as initial simulation time
% It0: distribution of H2BGFP intensities at initial time, used as prior (by default, log2 values; column vector) for random sampling
% lambda: average division rate (/week)
% indiv: number of independent simulations (number of simulated basal cells)
% tlag: refractory period between consecutive divisions (minimum cell-cycle period) (weeks)
% GamShape: 'Shape' parameter of the gamma-distributed cell-cycle period (=1 for the default exponential distribution)

%% Output:
% Idist: mxn matrix of H2BGFP intensities in n different basal cells across m time points

%% Example:
% rtime = [0 1 2 3]; %(weeks)
% It0 = 0; %presented as log2 value [original intensity = 1 or (2^0)]
% lambda = 2.9;
% indiv = 1000;
% tlag = 0.5/7;
% GamShape = 8;
% [Idist] = MonteCarloSimulator_BasalCell_H2BGFPdil(rtime,It0,lambda,indiv,tlag,GamShape);

%% Initial definition of parameters:
Idist = zeros(length(rtime),indiv);
GamScale = (1/(lambda)-tlag)/GamShape;

%% ITERATION FOR DIFFERENT INITIAL BASAL PROGENITOR CELLS
for aje = 1:indiv
    
    % Initial variables and cell attributes:
    time = rtime(1);
    It = 2.^(It0(randperm(size(It0,1),1),1)); % initial raw intensity (assumes prior It0 vector was given as a log2 value)(alternative: It = 2^(lognrnd(pd0.mu,pd0.sigma));)
    % Save the H2BGFP intensity at initial time point (t=0d): 
    Idist(1,aje) = It;
    
    % Asynchronous initial cond (defining a random value for the 'time-to-next-division' of initial cells at time 0)
    GamEstTime = [0:0.001:5]; % (wide time range for which to estimate the Gam-related integral)
    GamCumF = 1-gamcdf(GamEstTime,GamShape,GamScale);
    Prob_cutoff = tlag ./ (tlag + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
    myr1 = rand;
    if myr1 < Prob_cutoff
        tDiv0 = rand*tlag; % random sampling from a uniform dist.
    else
        tDiv0 = tlag + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.001; % random (multinomial) sampling from a gamma-related cum. dist.
    end
    
    % Calculate time to new division & Update H2BGFP content accordingly:
    time = time + tDiv0; %time = time + tlag+gamrnd(GamShape,GamScale); % update time
    for buc = 2:length(rtime)
        while time < rtime(buc)
            %It = It/2; % update intensity (assumption: non-random partitioning)
            %It = binornd(round(It/3),0.5)*3; % update intensity (assumption: random partitioning into daughter cells, amplified with dilution)
            It = It/2 + (rand-0.5)*sqrt(12)*0.06*It; % update intensity (assumption: random partitioning into daughter cells, with constant noise term)
            %constant noise term estimated from experiments c.a. 5% of It (1/sqrt(12) refers to the std of a uniform dist between [0,1])
            time = time + tlag+gamrnd(GamShape,GamScale);
        end
        % Save the H2BGFP intensity at subsequent time points (t=7d,12d,18d...): 
        Idist(buc,aje) = It;
    end
    
end
