function [nx_basal,nx_total,nx,ntime] = MonteCarloSimulator_SP_TotalCloneDynamics_Giroux2017(rtime,lambda,r,gamma,mu,indiv,tlag,GamShape)
%% Non-Markovian Monte Carlo simulator of Single Progenitor (SP) model dynamics - adapted for analysis of Giroux et al 2017 dataset
% Basal and total clone sizes are simulated over time starting from a single
% random labelled cell, which can be a basal or a suprabasal cell with a
% certain probability given by the initial distribution of GFP cells across
% layers at 24h post-induction (considered as t0).
% It allows to consider gamma-distributed cell cycle periods for progenitors,
% but also gamma-distributed periods for the shedding of suprabasal cells
% (i.e. the transition of differentiated cells through suprabasal layers 
% takes a minimum time).

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% lambda: average division rate (/week)
% r: probability of symmetric division leading to two progenitors, PP
% gamma: stratification rate (/week)
% mu: shedding rate (/week)
% indiv: number of independent simulations (number of simulated clones)
% tlag: refractory period between consecutive divisions (minimum cell-cycle period) or before shedding (in the case of suprabasal cells) (weeks)
% GamShape: 'Shape' parameter of the gamma-distributed cell-cycle period or, alternatively, time of residency in suprabasal layers (=1 for the default exponential distribution)

%% Output:
% nx_basal: mxn matrix of clone sizes (No. of basal cells per clone) over time (m clones x n time points)
% nx_total: mxn matrix of clone sizes (No. of total cells per clone) over time (m clones x n time points)
% nx: mxnxp matrix of clone sizes (No. of cells of each type per clone) over time (m clones x n time points x p cell types)
% ntime: horizontal vector of the n time points collected

%% Example:
% rtime = [1 4 26 52]; %(weeks)
% lambda = 2.9;
% r = 0.1;
% gamma = 5.4;
% mu = 1.5
% indiv = 1000;
% tlag = 0.5/7;
% GamShape = 8;
% [nx_basal,nx_total,nx,ntime] = MonteCarloSimulator_SP_TotalCloneDynamics_Giroux2017(rtime,lambda,r,gamma,mu,indiv,tlag,GamShape);

%% Default parameter values:
tic
% Initial error checks:
assert(nargin>5,'Single-Progenitor simulator function requires at least six inputs.')
assert((r>=0 & r<= 0.5),'Unvalid value for r. Should fall between 0-0.5')
if (nargin < 8)
    indiv=100;
    tlag = 0;
    GamShape = 1;
end
rho = gamma / (gamma + lambda); %proportion of progenitor basal cells (rho).

% Initial definition of parameters:
timelim=rtime(1,end); % time limit
ntime = rtime; % all collection time points
nx=zeros(indiv,length(rtime),3); % stores No. of cells of each type / clone
nx_basal=zeros(indiv,length(rtime)); % stores No. of basal cells / clone
nx_total=zeros(indiv,length(rtime)); % stores No. of total cells / clone
% GamScale (Gamma-dist scale param. that fits the observed average division rate - considering tlag)
GamScale = (1/lambda - tlag) ./ GamShape;
% GamScale (Gamma-dist scale param. that fits the observed average shedding rate - considering the tlag)
tlagShed = tlag;
GamShapeShed = GamShape;
GamScaleShed = (1/mu - tlagShed) ./ GamShapeShed;

%% ITERATION FOR DIFFERENT INDIVIDUAL CLONES
for it=1:indiv

    % Initial cell content of the clone (aprox. 90% of initial GFP cells are basal, the rest being suprabasal):
    %           P           D         SL
    x=mnrnd(1,[rho*0.9 (1-rho)*0.9 0.1]); % a single-cell clone with just one cell (can be progenitor / differentiating basal / suprabasal).

    % Initial variables & cell attributes:
    p=zeros(5,1);
    Id_Cell = 1; % Cell identifier within the clone
    tini = 0;
    tend = 0;
    
    % Asynchronous initial condition (defining a random value for the 'time-to-next-division' of initial cell at time 0)
    GamEstTime = [0:0.001:5]; % (wide time range for which to estimate the Gam-related integral)
    GamCumF = 1-gamcdf(GamEstTime,GamShape,GamScale);
    Prob_cutoff = tlag ./ (tlag + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
    myr1 = rand;
    if myr1 < Prob_cutoff
        tDiv0 = rand*tlag; % random sampling from a uniform dist.
    else
        tDiv0 = tlag + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.001; % random (multinomial) sampling from a gamma-related cum. dist.
    end 

    % Asynchronous initial cond (defining a random value for the 'time-to-next-shedding' of initial cells at time 0)
    GamEstTime = [0:0.001:5]; % (wide time range for which to estimate the Gam-related integral)
    GamCumF = 1-gamcdf(GamEstTime,GamShapeShed,GamScaleShed);
    Prob_cutoff = tlagShed ./ (tlagShed + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
    myr1 = rand;
    if myr1 < Prob_cutoff
        tShed0 = rand*tlagShed; % random sampling from a uniform dist.
    else
        tShed0 = tlagShed + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.001; % random (multinomial) sampling from a gamma-related cum. dist.
    end 

    % ITERATION FOR EACH SINGLE CELL WITHIN THE CLONE:
    while Id_Cell <= size(x,1)

        if tini(Id_Cell,1) <= timelim

            % Calculation of single-event probabilities:
            p(1) = lambda*r*x(Id_Cell,1); % P -> P + P
            p(2) = lambda*(1-2*r)*x(Id_Cell,1); % P -> P + D
            p(3) = lambda*r*x(Id_Cell,1); % P -> D + D
            p(4) = gamma*x(Id_Cell,2); % D -> c
            p(5) = mu*x(Id_Cell,3); % c -> lost
            % Calculation of total probability of event:
            pt=sum(p);

            % Calculate time to new event:
            if (x(Id_Cell,1)~=0) && (Id_Cell==1) % if founder progenitor
                tend(Id_Cell,1)=tini(Id_Cell,1)+tDiv0;
            elseif (x(Id_Cell,1)~=0) && (Id_Cell~=1) % any other progenitor
                tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShape,GamScale)+tlag;
            elseif (x(Id_Cell,3)~=0) && (Id_Cell==1) % if founder suprabasal
                tend(Id_Cell,1)=tini(Id_Cell,1)+tShed0;
            elseif (x(Id_Cell,3)~=0) && (Id_Cell~=1) % any other suprabasal
                tend(Id_Cell,1)=tini(Id_Cell,1)+gamrnd(GamShapeShed,GamScaleShed)+tlagShed;
            else % if differentiating basal cell (D)
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
                    x = [x; 0 0 1];
                    tini = [tini; tend(Id_Cell,1)];
%                 elseif (event==5) | loss is implicit in the time update
%                     x(3)=x(3)-1;
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

%% Sum all types of basal (P+D) and total (P+D+SL) cells to get basal-layer and total clone sizes:
nx_basal = nx(:,:,1)+nx(:,:,2);
nx_total = nx(:,:,1)+nx(:,:,2)+nx(:,:,3);
toc
