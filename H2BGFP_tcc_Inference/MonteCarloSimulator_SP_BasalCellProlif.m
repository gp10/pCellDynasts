function [NoDiv_BL] = MonteCarloSimulator_SP_BasalCellProlif(rtime,lambda,indiv,tlag,GamShape)
%% Non-Markovian Monte Carlo simulator of basal cell proliferation (number of division rounds over time)
% The number of division rounds experienced by a population of basal cells
% over time is simulated, as progenitor cells sucessively divide (only a
% one-branch pathway -the random walk from a given sibling- is considered
% upon each division).
% It allows to consider gamma-distributed cell cycle periods.

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% lambda: average division rate (/week)
% indiv: number of independent simulations (number of simulated basal cells)
% tlag: refractory period between consecutive divisions (minimum cell-cycle period) (weeks)
% GamShape: 'Shape' parameter of the gamma-distributed cell-cycle period (=1 for the default exponential distribution)

%% Output:
% NoDiv_BL: mxn matrix of No. of experienced rounds of divisions in n different basal cells across m time points

%% Example:
% rtime = [0 1 2 3]; %(weeks)
% lambda = 2.9;
% indiv = 1000;
% tlag = 0.5/7;
% GamShape = 8;
% [NoDiv_BL] = MonteCarloSimulator_SP_BasalCellProlif(rtime,lambda,indiv,tlag,GamShape);

%% Initial definition of parameters:
NoDiv_BL = zeros(length(rtime),indiv);
GamScale = (1/(lambda)-tlag)/GamShape;

%% ITERATION FOR DIFFERENT INITIAL BASAL PROGENITOR CELLS
for aje = 1:indiv
    
    % Initial variables and cell attributes:
    time = 0; % No of div at t=0d
    NoDiv = 0; % initial No. of cell division rounds experienced
    
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
    
    % Calculate time to new division & Update counter of cell division rounds accordingly:
    time = time + tDiv0; %time = time + tlag+gamrnd(GamShape,GamScale); % update time
    for buc = 1:length(rtime)
        while time < rtime(buc)
            %It = It/2; % update intensity (assumption: non-random partitioning)
            %It = binornd(round(It/3),0.5)*3; % update intensity (assumption: random partitioning into daughter cells, amplified with dilution)
            %It = It/2 + (rand-0.5)*sqrt(12)*0.06*It; % update intensity (assumption: random partitioning into daughter cells, with constant noise term)
            %constant noise term estimated from experiments c.a. 5% of It (1/sqrt(12) refers to the std of a uniform dist between [0,1])
            NoDiv = NoDiv + 1;
            time = time + tlag+gamrnd(GamShape,GamScale);
        end
        % Save the No. of cell division rounds experienced by each time point (t=3d,6d,9d...): 
        NoDiv_BL(buc,aje) = NoDiv;
    end
    
end
