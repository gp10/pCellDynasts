function [nx_basal,ntime] = MonteCarloSimulator_SP_BasalCloneDynamics_Fullgrabe2015(rtime,iniCloneSizes,lambda,r,gamma,indiv,tlag,GamShape)
%% Non-Markovian Monte Carlo simulator of Single Progenitor (SP) model dynamics - adapted for Füllgrabe et al 2015 data analysis
% Basal clone sizes are simulated over time starting from preformed clones
% consisting of multiple (progenitor and/or differentiating) labelled cells.
% It allows to consider gamma-distributed cell cycle periods.

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% iniCloneSizes: vector of [NClones,1] shape containing initial experimental clone sizes
% lambda: average division rate (/week)
% r: probability of symmetric division leading to two progenitors, PP
% gamma: stratification rate (/week)
% indiv: number of independent simulations (number of simulated clones)
% tlag: refractory period between consecutive divisions (minimum cell-cycle period) (weeks)
% GamShape: 'Shape' parameter of the gamma-distributed cell-cycle period (=1 for the default exponential distribution)

%% Output:
% nx_basal: mxn matrix of clone sizes (No. of basal cells per clone) over time (m clones x n time points)
% ntime: horizontal vector of the n time points collected

%% Example:
% rtime = [1 4 26 52]; %(weeks)
% lambda = 2.9;
% r = 0.1;
% gamma = 5.4;
% indiv = 1000;
% tlag = 0.5/7;
% GamShape = 8;
% [nx_basal,ntime] = MonteCarloSimulator_SP_BasalCloneDynamics_Fullgrabe2015(rtime,iniCloneSizes,lambda,r,gamma,indiv,tlag,GamShape);

%% Default parameter values:
tic
% Initial error checks:
assert(nargin>3,'Single-Progenitor simulator function requires at least four inputs.')
assert((r>=0 & r<= 0.5),'Unvalid value for r. Should fall between 0-0.5')
if (nargin < 7)
    indiv=100;
    tlag = 0;
    GamShape = 1;
end
rho = gamma / (gamma + lambda); %proportion of progenitor basal cells (rho).

% Initial definition of parameters:
timelim=rtime(1,end); % time limit
ntime = rtime; % all collection time points
nx=zeros(indiv,length(rtime),2); % stores No. of cells of each type / clone
nx_basal=zeros(indiv,length(rtime)); % stores No. of basal cells / clone
% GamScale (Gamma-dist scale param. that fits the observed average division rate - considering tlag)
GamScale = (1/lambda - tlag) ./ GamShape;

%% ITERATION FOR DIFFERENT INDIVIDUAL CLONES
for it=1:indiv

    % Initial cell content of the clone:
    %  P  D
    %x=[1  0]; % a single-cell clone with just one progenitor cell (P).
    [x,tini,tend,iniCloneSize] = sampling_Initial_Clone_Fullgrabe2015(iniCloneSizes,rho); % multicellular clone of random P,D cell composition at the beginning

    % Initial variables & cell attributes:
    p=zeros(4,1);
    Id_Cell = 1; % Cell identifier within the clone
    %tini = 0;
    %tend = 0;
    
    % Asynchronous initial condition (defining a random value for the 'time-to-next-division' of initial cells at time 0)
    GamEstTime = [0:0.001:5]; % (wide time range for which to estimate the Gam-related integral)
    GamCumF = 1-gamcdf(GamEstTime,GamShape,GamScale);
    Prob_cutoff = tlag ./ (tlag + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
    % code continued inside the loop... 

    % ITERATION FOR EACH SINGLE CELL WITHIN THE CLONE:
    while Id_Cell <= size(x,1)

        if tini(Id_Cell,1) <= timelim

            % Calculation of single-event probabilities:
            p(1) = lambda*r*x(Id_Cell,1); % P -> P + P
            p(2) = lambda*(1-2*r)*x(Id_Cell,1); % P -> P + D
            p(3) = lambda*r*x(Id_Cell,1); % P -> D + D
            p(4) = gamma*x(Id_Cell,2); % D -> lost
            % Calculation of total probability of event:
            pt=sum(p);

            % Calculate time to new event:
            if (x(Id_Cell,1)~=0) && (Id_Cell<=iniCloneSize) % if preexistent progenitor
                myr1 = rand;
                if myr1 < Prob_cutoff
                    tDiv0 = rand*tlag; % random sampling from a uniform dist.
                else
                    tDiv0 = tlag + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.001; % random (multinomial) sampling from a gamma-related cum. dist.
                end
                tend(Id_Cell,1)=tini(Id_Cell,1)+tDiv0;
            elseif (x(Id_Cell,1)~=0) && (Id_Cell>iniCloneSize) % any other progenitor
                tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShape,GamScale)+tlag;
            else % if differentiating cell (D)
                tau=-(1./pt)*log(rand);
                tend(Id_Cell,1)=tini(Id_Cell,1)+tau;
            end

            if tend(Id_Cell,1) < timelim
                % Event selection:
                event = find(cumsum(p)>(rand*pt),1);
                if (event==1)
                    x = [x; 1 0; 1 0];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                elseif (event==2)
                    x = [x; 1 0; 0 1];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                elseif (event==3)
                    x = [x; 0 1; 0 1];
                    tini = [tini; tend(Id_Cell,1); tend(Id_Cell,1)];
                %elseif (event==4) | loss is implicit in the time update
                end
            end

        end
        Id_Cell = Id_Cell + 1;

    end

    % Save the populations of cells at indicated time points:
    for bas = 1:size(rtime,2)
        nx(it,bas,1) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,1)==1 ),1);
        nx(it,bas,2) = size(find( (tini <= rtime(1,bas)) & (rtime(1,bas) < tend) & x(:,2)==1 ),1);
    end

end

%% Sum both types of basal cells (P+D) to get basal-layer clone sizes:
nx_basal = nx(:,:,1)+nx(:,:,2);
toc
