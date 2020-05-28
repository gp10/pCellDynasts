function [NoDiv_SL] = MonteCarloSimulator_SP_SuprabasalCellDiv(rtime,indiv,rho,lambda,tlagDiv,GamShapeDiv,mu,tlagShed,GamShapeShed)
%% Non-Markovian Monte Carlo simulator of suprabasal cell replacement by basal proliferation (number of division rounds over time)
% The number of division rounds experienced by suprabasal cells is
% simulated over time, as progenitor cells sucessively divide in the basal
% layer (only a one-branch pathway -the random walk from a given sibling-
% is considered upon each division) and a daughter stratifies into the
% spinous (suprabasal) layer.
% It allows to consider gamma-distributed cell cycle and stratification/sheeding periods.

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing desired time points (expressed in weeks)
% indiv: number of independent simulations (number of simulated spinous suprabasal cells)
% rho: proportion of basal cell being dividing progenitors
% lambda: average division rate (/week)
% tlagDiv: refractory period between consecutive divisions (minimum cell-cycle period) (weeks)
% GamShapeDiv: 'Shape' parameter of the gamma-distributed cell-cycle period (=1 for the default exponential distribution)
% mu: average shedding rate from spinous suprabasal layer (/week)
% tlagShed: refractory period before shedding of spinous suprabasal cells (minimum time for transitioning through suprbasal layers) (weeks)
% GamShapeShed: 'Shape' parameter of the gamma-distributed waiting time period for shedding (=1 for the default exponential distribution)

%% Output:
% NoDiv_SL: mxn matrix of No. of experienced rounds of divisions in n different spinous (suprabasal) cells across m time points

%% Example:
% rtime = [0 1 2 3]; %(weeks)
% indiv = 1000;
% rho = 0.5;
% lambda = 2.9;
% tlagDiv = 0.5/7;
% GamShapeDiv = 8;
% mu = 2.9;
% tlagShed = 0.5/7;
% GamShapeShed = 8;
% [NoDiv_SL] = MonteCarloSimulator_SP_SuprabasalCellDiv(rtime,indiv,rho,lambda,tlagDiv,GamShapeDiv,mu,tlagShed,GamShapeShed);

%% Initial definition of parameters:
NoDiv_SL = zeros(length(rtime),indiv);
GamScaleDiv = (1/(lambda)-tlagDiv)/GamShapeDiv;
GamScaleShed = (1/(mu)-tlagShed)/GamShapeShed;

%% ITERATION FOR DIFFERENT SPINOUS (SUPRABASAL) CELLS
for aje = 1:indiv
    
    for buc = 1:length(rtime)
    
        % Initial variables and cell attributes:
        time = 0; % No of div at t=0d
        NoDiv = 0; % initial No. of cell division rounds experienced
        layer = 1; % [0=BL | 1=SL] (i.e. we start with a single SL cell)
        timeShed = 0; % time for shedding (default value)
        
        while time < rtime(buc)

            % Simulate BL cell:
            while (layer==0)
                % progenitor cell (A):
                if (rand <= rho)
                    if (time == 0)
                        % Asynchronous initial cond (defining a random value for the 'time-to-next-division' of initial cells at time 0)
                        GamEstTime = [0:0.001:5]; % (wide time range for which to estimate the Gam-related integral)
                        GamCumF = 1-gamcdf(GamEstTime,GamShapeDiv,GamScaleDiv);
                        Prob_cutoff = tlagDiv ./ (tlagDiv + 0.001.*trapz(GamCumF)); % defined by the density of each of the 2 cdfs
                        myr1 = rand;
                        if myr1 < Prob_cutoff
                            tDiv0 = rand*tlagDiv; % random sampling from a uniform dist.
                        else
                            tDiv0 = tlagDiv + GamEstTime(1,find(mnrnd(1,GamCumF./sum(GamCumF))==1)) + rand*0.001; % random (multinomial) sampling from a gamma-related cum. dist.
                        end
                        time = time + tDiv0; %time = time + tlag+gamrnd(GamShape,GamScale); % update time
                    else
                        time = time + tlagDiv+gamrnd(GamShapeDiv,GamScaleDiv);
                    end
                    if (time >= timeShed) % If cell A persists until timeShed, it is forced to migrate to SL (time is reset to simulate that SL from then onwards)
                        time = timeShed;
                        layer = 1;
                    else % If cell A divides before timeShed, NoDiv increases and a random BL cell is simulated from then onwards
                        layer = 0;
                        NoDiv = NoDiv + 1;
                    end
                % differentiating cell (B):
                else
                    if (time == 0)
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
                        time = time + tShed0; %time = time + tlag+gamrnd(GamShape,GamScale); % update time
                    else
                        time = time + tlagShed+gamrnd(GamShapeShed,GamScaleShed);
                    end
                    if (time >= timeShed) % If cell B persists until timeShed, it is forced to migrate to SL (time is reset to simulate that SL from then onwards)
                        time = timeShed;
                        layer = 1;
                    else % If cell B leaves BL before timeShed, it is not a good replacement, i.e. it is ignored and we restart with a random BL cell from the beginning
                        layer = 0;
                        time = 0;
                        NoDiv = 0;
                    end
                end
            end

            % Simulate SL cell:
            if (layer == 1)
                if (time == 0)
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
                    time = time + tShed0; %time = time + tlag+gamrnd(GamShape,GamScale); % update time
                    timeShed = time;
                else
                    time = time + tlagShed+gamrnd(GamShapeShed,GamScaleShed);
                    timeShed = time;
                end
                if time < rtime(buc) % If cell leaves SL before acquisition time, gap is to be replaced by BL (to be simulated from the beginning)
                    layer = 0;
                    time = 0;
                    NoDiv = 0;
                end % If cell persists until acquisition time, it is recorded
            end
                        
        end
        % Save the No. of cell division rounds experienced by each time point (t=3d,6d,9d...):
        NoDiv_SL(buc,aje) = NoDiv;
    
    end
    
end
