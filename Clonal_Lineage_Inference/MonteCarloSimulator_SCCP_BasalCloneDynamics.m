function [nx_basal,ntime] = MonteCarloSimulator_SCCP_BasalCloneDynamics(rtime,labelSfrac,lambdaS,rS,lambdaP,rP,DeltaP,gamma,indiv,tlagS,GamShapeS,tlagP,GamShapeP)
%% Non-Markovian Monte Carlo simulator of Stem Cell-Committed Progenitor (SC-CP) model dynamics
% Basal clone sizes are simulated over time starting from a single labelled
% proliferating cell. It allows to consider gamma-distributed cell cycle periods.

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% labelSfrac: fraction of initially labelled cells being stem cells
% lambdaS: average division rate of S-population (/week)
% rS: probability of symmetric division of S leading to two stem cells, SS
% lambdaP: average division rate of P-population (/week)
% rP: probability of symmetric division of P leading to two progenitors, PP (in the absence of fate imbalance)
% DeltaP: fate imbalance in P division outcome
% gamma: stratification rate (/week)
% indiv: number of independent simulations (number of simulated clones)
% tlagS: refractory period between consecutive S divisions (minimum cell-cycle period for S cells) (weeks)
% GamShapeS: 'Shape' parameter of the gamma-distributed cell-cycle period for S cells (=1 for the default exponential distribution)
% tlagP: refractory period between consecutive P divisions (minimum cell-cycle period for P cells) (weeks)
% GamShapeP: 'Shape' parameter of the gamma-distributed cell-cycle period for P cells (=1 for the default exponential distribution)

%% Output:
% nx_basal: mxn matrix of clone sizes (No. of basal cells per clone) over time (m clones x n time points)
% ntime: horizontal vector of the n time points collected

%% Example:
% rtime = [1 4 26 52]; %(weeks)
% labelSfrac = 0.65;
% lambdaS = 0.45;
% rS = 0.03;
% lambdaP = 1.7;
% rP = 0.19;
% DeltaP = 0.02;
% gamma = 118.8;
% indiv = 1000;
% tlagS = 0/7;
% GamShapeS = 1;
% tlagP = 0/7;
% GamShapeP = 1;
% [nx_basal,ntime] = MonteCarloSimulator_SCCP_BasalCloneDynamics(rtime,labelSfrac,lambdaS,rS,lambdaP,rP,DeltaP,gamma,indiv,tlagS,GamShapeS,tlagP,GamShapeP);

%% Default parameter values:
tic
% Initial error checks:
assert(nargin>7,'SC-CP model simulator function requires at least eight inputs.')
assert((rS>=0 & rS<= 0.5),'Unvalid value for rS. Should fall between 0-0.5')
assert((rP>=0 & rP<= 0.5),'Unvalid value for rP. Should fall between 0-0.5')
if (nargin < 9)
    indiv=100;
    tlagS = 0;
    GamShapeS = 1;
    tlagP = 0;
    GamShapeP = 1;
end

% Initial definition of parameters:
timelim=rtime(1,end); % time limit
ntime = rtime; % all collection time points
nx=zeros(indiv,length(rtime),3); % stores No. of cells of each type / clone
nx_basal=zeros(indiv,length(rtime)); % stores No. of basal cells / clone
% GamScale (Gamma-dist scale param. that fits the observed average division rate - considering tlag)
GamScaleP = (1/lambdaP - tlagP) ./ GamShapeP;
GamScaleS = (1/lambdaS - tlagS) ./ GamShapeS;

%% ITERATION FOR DIFFERENT INDIVIDUAL CLONES
for it=1:indiv

    % Initial cell content of the clone:
    %                 S         P       D
    x= mnrnd(1,[labelSfrac 1-labelSfrac 0]); % a single-cell clone with with either one stem cell or one progenitor (proportion given by labelSfrac).

    % Initial variables & cell attributes:
    p=zeros(7,1);
    Id_Cell = 1; % Cell identifier within the clone
    tini = 0;
    tend = 0;
    
    % Asynchronous initial cond (defining a random value for the 'time-to-next-division' of initial S cells at time 0)
    GamEstTime = [0:0.01:20]; % (wide time range for which to estimate the Gam-related integral)
    GamCumF = 1-gamcdf(GamEstTime,GamShapeS,GamScaleS);
    Prob_cutoff = tlagS ./ (tlagS + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
    myr1 = rand;
    if myr1 < Prob_cutoff
        tDiv0_S = rand*tlagS; % random sampling from a uniform dist.
    else
        tDiv0_S = tlagS + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.01; % random (multinomial) sampling from a gamma-related cum. dist.
    end

    % Asynchronous initial cond (defining a random value for the 'time-to-next-division' of initial P cells at time 0)
    GamEstTime = [0:0.001:5]; % (wide time range for which to estimate the Gam-related integral)
    GamCumF = 1-gamcdf(GamEstTime,GamShapeP,GamScaleP);
    Prob_cutoff = tlagP ./ (tlagP + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
    myr1 = rand;
    if myr1 < Prob_cutoff
        tDiv0_P = rand*tlagP; % random sampling from a uniform dist.
    else
        tDiv0_P = tlagP + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.001; % random (multinomial) sampling from a gamma-related cum. dist.
    end

    % ITERATION FOR EACH SINGLE CELL WITHIN THE CLONE:
    while Id_Cell <= size(x,1)

        if tini(Id_Cell,1) <= timelim

            % Calculation of single-event probabilities:
            p(1) = lambdaS*rS*x(Id_Cell,1); % S -> S + S
            p(2) = lambdaS*(1-2*rS)*x(Id_Cell,1); % S -> S + P
            p(3) = lambdaS*rS*x(Id_Cell,1); % S -> P + P
            p(4) = lambdaP*(1-DeltaP)*rP*x(Id_Cell,2); % P -> P + P
            p(5) = lambdaP*(1-2*rP)*x(Id_Cell,2); % P -> P + D
            p(6) = lambdaP*(1+DeltaP)*rP*x(Id_Cell,2); % P -> D + D
            p(7) = gamma*x(Id_Cell,3); % D -> lost
            % Calculation of total probability of event:
            pt=sum(p);

            % Calculate time to new event:
            if (x(Id_Cell,1)~=0) && (Id_Cell==1) % if founder stem cell
                tend(Id_Cell,1)=tini(Id_Cell,1)+tDiv0_S;
            elseif (x(Id_Cell,1)~=0) && (Id_Cell~=1) % any other stem cell
                tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShapeS,GamScaleS)+tlagS;
            elseif (x(Id_Cell,2)~=0) && (Id_Cell==1) % if founder progenitor
                tend(Id_Cell,1)=tini(Id_Cell,1)+tDiv0_P;
            elseif (x(Id_Cell,2)~=0) && (Id_Cell~=1) % any other progenitor
                tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShapeP,GamScaleP)+tlagP;
            else % if differentiating cell (D)
                tau=-(1./pt)*log(rand);
                tend(Id_Cell,1)=tini(Id_Cell,1)+tau;
            end

            if tend(Id_Cell,1) < timelim
                % Event selection:
                event = find(cumsum(p)>(rand*pt),1);
                if (event==1)
                    x = [x; 1 0 0; 1 0 0];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                elseif (event==2)
                    x = [x; 1 0 0; 0 1 0];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                elseif (event==3)
                    x = [x; 0 1 0; 0 1 0];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                elseif (event==4)
                    x = [x; 0 1 0; 0 1 0];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                elseif (event==5)
                    x = [x; 0 1 0; 0 0 1];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                elseif (event==6)
                    x = [x; 0 0 1; 0 0 1];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
%                 elseif (event==7) | loss is implicit in the time update
                end
            end

        end
        Id_Cell = Id_Cell + 1;

    end

    % Save the populations of cells at indicated time points:
    for bas = 1:size(rtime,2)
        nx(it,bas,1) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,1)==1 ),1);
        nx(it,bas,2) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,2)==1 ),1);
        nx(it,bas,3) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,3)==1 ),1);
    end

end

%% Sum all types of basal cells (S+P+D) to get basal-layer clone sizes:
nx_basal = nx(:,:,1)+nx(:,:,2)+nx(:,:,3);
toc
