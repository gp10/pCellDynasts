function [nx_basal,nx_total,ntime] = MonteCarloSimulator_2xSC_TotalCloneDynamics(rtime,labelS1frac,lambdaS1,lambdaS2,uS1,gammaS1,gammaS2,mu,indiv,tlagDivS1,GamShapeDivS1,tlagStrS1,GamShapeStrS1,tlagDivS2,GamShapeDivS2,tlagStrS2,GamShapeStrS2,tlagShed,GamShapeShed)
%% Non-Markovian Monte Carlo simulator of 2 Independent Stem Cell (2xSC) model dynamics
% Basal clone sizes are simulated over time starting from a single labelled
% proliferating cell. It allows to consider gamma-distributed cell cycle periods.

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% labelS1frac: fraction of initially labelled cells being type-1 stem cells, S1
% lambdaS1: average division rate of S1 population (/week)
% lambdaS2: average division rate of S2 population (/week)
% uS1: probability that S1 stratification occurs uncoupled from division
% gammaS1: division-uncoupled stratification rate for S1 (/week)
% gammaS2: division-uncoupled stratification rate for S2 (/week)
% mu: shedding rate (/week) or spinous-to-granular layer stratification rate
% indiv: number of independent simulations (number of simulated clones)
% tlagDivS1: refractory period between consecutive S1 divisions (minimum cell-cycle period for S1 cells) (weeks)
% GamShapeDivS1: 'Shape' parameter of the gamma-distributed cell-cycle period for dividing S1 cells (=1 for the default exponential distribution)
% tlagStrS1: refractory period for the time of a newborn S1 to stratify (independently of a division event) (weeks)
% GamShapeStrS1: 'Shape' parameter of the gamma-distributed time for division-uncoupled stratification of S1 cells (=1 for the default exponential distribution)
% tlagDivS2: refractory period between consecutive S2 divisions (minimum cell-cycle period for S2 cells) (weeks)
% GamShapeDivS2: 'Shape' parameter of the gamma-distributed cell-cycle period for dividing S2 cells (=1 for the default exponential distribution)
% tlagStrS2: refractory period for the time of a newborn S2 to stratify (weeks)
% GamShapeStrS2: 'Shape' parameter of the gamma-distributed time for stratification of S2 cells (=1 for the default exponential distribution)
% tlagShed: refractory period for the time of a new suprabasal (spinous) cell to shed (stratify to granular layer) (weeks)
% GamShapeShed: 'Shape' parameter of the gamma-distributed time for shedding (or spinous-to-granular transition) of suprabasal cells (=1 for the default exponential distribution)

%% Output:
% nx_basal: mxn matrix of clone sizes (No. of basal cells per clone) over time (m clones x n time points)
% nx_total: mxn matrix of clone sizes (No. of total cells per clone) over time (m clones x n time points)
% ntime: horizontal vector of the n time points collected

%% Example:
% rtime = [1 4 26 52]; %(weeks)
% labelS1frac = 0.70; % arbitrary
% lambdaS1 = 0.51*7; %/week
% lambdaS2 = 0.19*7; %/week
% uS1 = 0;
% gammaS1 = lambdaS1;
% gammaS2 = lambdaS2;
% mu = 3.4596; % /week
% indiv = 1000;
% tlagDivS1 = 0;
% GamShapeDivS1 = 2;
% tlagStrS1 = 0;
% GamShapeStrS1 = 2;
% tlagDivS2 = 0;
% GamShapeDivS2 = 1;
% tlagStrS2 = 0;
% GamShapeStrS2 = 1;
% tlagShed = 0;
% GamShapeShed = 2;
%[nx_basal,nx_total,ntime] = MonteCarloSimulator_2xSC_TotalCloneDynamics(rtime,labelS1frac,lambdaS1,lambdaS2,uS1,gammaS1,gammaS2,mu,indiv,tlagDivS1,GamShapeDivS1,tlagStrS1,GamShapeStrS1,tlagDivS2,GamShapeDivS2,tlagStrS2,GamShapeStrS2,tlagShed,GamShapeShed);

%% Default parameter values:
tic
% Initial error checks:
assert(nargin>7,'2xSC model simulator function requires at least eight inputs.')
assert((uS1>=0 & uS1<= 1),'Unvalid value for uS1. Should fall between 0-1')
if (nargin < 19)
    indiv=100;
    tlagDivS1 = 0;
    GamShapeDivS1 = 1;
    tlagStrS1 = 0;
    GamShapeStrS1 = 1;
    tlagDivS2 = 0;
    GamShapeDivS2 = 1;
    tlagStrS2 = 0;
    GamShapeStrS2 = 1;
    tlagShed = 0;
    GamShapeShed = 1;
end

% Initial definition of parameters:
timelim=rtime(1,end); % time limit
ntime = rtime; % all collection time points
nx=zeros(indiv,length(rtime),3); % stores No. of cells of each type / clone
nx_basal=zeros(indiv,length(rtime)); % stores No. of basal cells / clone
nx_total=zeros(indiv,length(rtime)); % stores No. of total cells / clone
% GamScale (Gamma-dist scale param. that fits the observed average division/stratification/shedding rate - considering the respective tlag)
GamScaleDivS1 = (1/lambdaS1 - tlagDivS1) ./ GamShapeDivS1;
GamScaleStrS1 = (1/gammaS1 - tlagStrS1) ./ GamShapeStrS1;
GamScaleDivS2 = (1/lambdaS2 - tlagDivS2) ./ GamShapeDivS2;
GamScaleStrS2 = (1/gammaS2 - tlagStrS2) ./ GamShapeStrS2;
GamScaleShed = (1/mu - tlagShed) ./ GamShapeShed;

%% ITERATION FOR DIFFERENT INDIVIDUAL CLONES
for it=1:indiv

    % Initial cell content of the clone:
    %               S1          S2        SL
    x= mnrnd(1,[labelS1frac 1-labelS1frac 0]); % a single-cell clone with either one stem cell typeI or one stem cell typeII (proportion given by labelS1frac).

    % Initial variables & cell attributes:
    p=zeros(6,1);
    Id_Cell = 1; % Cell identifier within the clone
    tini = 0;
    tend = 0;

    % Asynchronous initial cond (defining a random value for the 'time-to-next-division' of initial S1 cells at time 0)
    GamEstTime = [0:0.01:20]; % (wide time range for which to estimate the Gam-related integral)
    GamCumF = 1-gamcdf(GamEstTime,GamShapeDivS1,GamScaleDivS1);
    Prob_cutoff = tlagDivS1 ./ (tlagDivS1 + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
    myr1 = rand;
    if myr1 < Prob_cutoff
        tDiv0_S1 = rand*tlagDivS1; % random sampling from a uniform dist.
    else
        tDiv0_S1 = tlagDivS1 + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.01; % random (multinomial) sampling from a gamma-related cum. dist.
    end

    % Asynchronous initial cond (defining a random value for the 'time-to-next-stratification' of initial S1 cells at time 0)
    GamEstTime = [0:0.01:20]; % (wide time range for which to estimate the Gam-related integral)
    GamCumF = 1-gamcdf(GamEstTime,GamShapeStrS1,GamScaleStrS1);
    Prob_cutoff = tlagStrS1 ./ (tlagStrS1 + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
    myr1 = rand;
    if myr1 < Prob_cutoff
        tStr0_S1 = rand*tlagStrS1; % random sampling from a uniform dist.
    else
        tStr0_S1 = tlagStrS1 + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.01; % random (multinomial) sampling from a gamma-related cum. dist.
    end

    % Asynchronous initial cond (defining a random value for the 'time-to-next-division' of initial S2 cells at time 0)
    GamEstTime = [0:0.001:5]; % (wide time range for which to estimate the Gam-related integral)
    GamCumF = 1-gamcdf(GamEstTime,GamShapeDivS2,GamScaleDivS2);
    Prob_cutoff = tlagDivS2 ./ (tlagDivS2 + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
    myr1 = rand;
    if myr1 < Prob_cutoff
        tDiv0_S2 = rand*tlagDivS2; % random sampling from a uniform dist.
    else
        tDiv0_S2 = tlagDivS2 + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.001; % random (multinomial) sampling from a gamma-related cum. dist.
    end

    % Asynchronous initial cond (defining a random value for the 'time-to-next-stratification' of initial S2 cells at time 0)
    GamEstTime = [0:0.001:5]; % (wide time range for which to estimate the Gam-related integral)
    GamCumF = 1-gamcdf(GamEstTime,GamShapeStrS2,GamScaleStrS2);
    Prob_cutoff = tlagStrS2 ./ (tlagStrS2 + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
    myr1 = rand;
    if myr1 < Prob_cutoff
        tStr0_S2 = rand*tlagStrS2; % random sampling from a uniform dist.
    else
        tStr0_S2 = tlagStrS2 + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.001; % random (multinomial) sampling from a gamma-related cum. dist.
    end

    % ITERATION FOR EACH SINGLE CELL WITHIN THE CLONE:
    while Id_Cell <= size(x,1)

        if tini(Id_Cell,1) <= timelim

            % Calculation of single-event probabilities:
            p(1) = lambdaS1*uS1/(1+uS1)*x(Id_Cell,1); % S1 -> S1 + S1
            p(2) = gammaS1*uS1/(1+uS1)*x(Id_Cell,1); % S1 -> c
            p(3) = lambdaS1*(1-uS1)/(1+uS1)*x(Id_Cell,1); % S1 -> S1 + c
            p(4) = lambdaS2*0.5*x(Id_Cell,2); % S2 -> S2 + S2
            p(5) = gammaS2*0.5*x(Id_Cell,2); % S2 -> c
            p(6) = mu*x(Id_Cell,3); % c -> lost
            % Calculation of total probability of event:
            pt=sum(p);

            % Calculate time to new event:
            % (specified below)

            % Event selection:
            event = find(cumsum(p)>(rand*pt),1);
            if (event==1)
                if Id_Cell==1 % if founder type-1 stem cell
                    tend(Id_Cell,1)=tini(Id_Cell,1)+tDiv0_S1;
                else % any other type-1 stem cell
                    tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShapeDivS1,GamScaleDivS1)+tlagDivS1;
                end
                if tend(Id_Cell,1) < timelim
                    x = [x; 1 0 0; 1 0 0];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                end
            elseif (event==2)
                if Id_Cell==1 % if founder type-1 stem cell
                    tend(Id_Cell,1)=tini(Id_Cell,1)+tStr0_S1;
                else % any other type-1 stem cell
                    tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShapeStrS1,GamScaleStrS1)+tlagStrS1;
                end
                if tend(Id_Cell,1) < timelim
                    x = [x; 0 0 1];
                    tini = [tini; tend(Id_Cell,1)];
                end
            elseif (event==3) % if founder type-1 stem cell
                if Id_Cell==1
                    tend(Id_Cell,1)=tini(Id_Cell,1)+tDiv0_S1;
                else % any other type-1 stem cell
                    tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShapeDivS1,GamScaleDivS1)+tlagDivS1;
                end
                if tend(Id_Cell,1) < timelim
                    x = [x; 1 0 0; 0 0 1];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                end
            elseif (event==4)
                if Id_Cell==1 % if founder type-2 stem cell
                    tend(Id_Cell,1)=tini(Id_Cell,1)+tDiv0_S2;
                else % any other type-2 stem cell
                    tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShapeDivS2,GamScaleDivS2)+tlagDivS2;
                end
                if tend(Id_Cell,1) < timelim
                    x = [x; 0 1 0; 0 1 0];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                end
            elseif (event==5)
                if Id_Cell==1 % if founder type-2 stem cell
                    tend(Id_Cell,1)=tini(Id_Cell,1)+tStr0_S2;
                else % any other type-2 stem cell
                    tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShapeStrS2,GamScaleStrS2)+tlagStrS2;
                end
                if tend(Id_Cell,1) < timelim
                    x = [x; 0 0 1];
                    tini = [tini; tend(Id_Cell,1)];
                end
            elseif (event==6) % loss is implicit in the time update, but not its timing!
                tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShapeShed,GamScaleShed)+tlagShed;
            end
            
        end
        Id_Cell = Id_Cell + 1;

    end

    % Save the populations of cells at certain time points:
    for bas = 1:size(rtime,2)
        nx(it,bas,1) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,1)==1 ),1);
        nx(it,bas,2) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,2)==1 ),1);
        nx(it,bas,3) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,3)==1 ),1);
    end

end

%% Sum all types of basal (S1+S2) and total (S1+S2+SL) cells to get basal-layer and total clone sizes:
nx_basal = nx(:,:,1)+nx(:,:,2);
nx_total = nx(:,:,1)+nx(:,:,2)+nx(:,:,3);
toc
