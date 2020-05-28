%% SCRIPT USED FOR tcc PARAMETER INFERENCE AND FOR SP MODEL GOODNESS-OF-FIT ON DISTRIBUTIONS OF CELL PROLIFERATION FROM MASCRE ET AL (2012):
% The average division rate is calculated and an Approximated Bayesian
% Computation (ABC) algorithm can be run to infer the distribution of
% individual progenitor tcc that best suits the experimental H2BGFP-
% deconvoluted histograms of the No. of division rounds shown for basal
% cells by Mascre et al (2012).

%% SELECTION OF GENERAL SETTINGS:
runABC = 0; % ( 1=runs ABC algorithm to compute likeliest SP cell-cycle period distributions matching cell division profiles | 0=skips ABC calculations and loads already preset optimal cell-cycle period distributions referred to in the text )
set_TolThres = 1; % ( 1=Run a test to set an adequate ABC tolerance threshold for a given acceptance rate | 0=Use preset adequate tolerance threshold levels )

%% LOAD H2BGFP-DECONVOLUTED DATA ON FREQUENCIES OF BASAL CELL DIVISION NUMBER AT DIFFERENT TIMES FROM Mascre et al (2012):
load ./Datasets/H2BGFPdil_Tail_Mascre2012_dataset.mat

%% PLOT REPORTED HISTOGRAMS OF No. OF CELL DIVISION ROUNDS:
figure(1)
% preset parameters (as in the original paper):
xtickpos = 1:10;
xlabelType = {'0','1','2','3','4','5','6','7','8','9'};
myylim = [0 0.8; 0 0.8; 0 0.8; 0 0.4; 0 0.4; 0 0.4];
rtime = ntime;

% sketch colors (just for legend-display purposes):
ColPalette = 'br'; StylePalette = {'-','-'};
subplot(3,3,1); hold on; for palb = 1:2; plot([0 0],[0 0],'Color',ColPalette(palb),'LineStyle',StylePalette{1,palb}); end
legend('SC-CP model fit (Mascré et al)','SP model fit (optim param.)')

for aja = 1:length(ntime) %different time points
    % BL
    subplot(3,3,aja)
    hold on
    myerr(:,:,1) = abs(FreqDiv_BL_mean{1,aja}-FreqDiv_BL_dn{1,aja});
    myerr(:,:,2) = abs(FreqDiv_BL_up{1,aja}-FreqDiv_BL_mean{1,aja});
    barwitherr(myerr,FreqDiv_BL_mean{1,aja},'FaceColor','g'); ylim(myylim(aja,:))
    xlim([xtickpos(1)-0.7 xtickpos(end)+0.7]); set(gca,'XTick',xtickpos); set(gca,'XTickLabel',xlabelType);
    xlabel('No. cell divisions'); ylabel('Frequency')
end

%% PLOT MASCRE ET AL (2012) SC-CP MODEL FIT:
for aja = 1:length(ntime) % different time points
    subplot(3,3,aja)
    hold on
    plot(FreqDiv_BL_SCCPmodel{1,aja},'b');
end

%% ESTIMATE AVERAGE DIVISION RATE
% Retrieve average number of divisions over time:
avg_No_div = [];
for aja = 1:length(ntime)
    avg_No_div(1,aja) = sum( FreqDiv_BL_mean{1,aja}.* (xtickpos-1)' ) ./ sum(FreqDiv_BL_mean{1,aja});
end
% Linear fit on the average number of divisions over time:
myfun = fittype('p1*x + p2','dependent',{'y'},'independent',{'x'},'coefficients',{'p1', 'p2'});
[LineFit_all,LineFit_all_stat] = fit(ntime',avg_No_div',myfun,'Startpoint',[0 0],'Lower',[-Inf 0],'Upper',[Inf 0]);
Lambda_avg = abs(LineFit_all.p1);
disp(sprintf('Avg. Division rate: lambda = %.2f/week',Lambda_avg));

% Plot best fit for the average division rate:
figure()
hold on
plot(ntime', avg_No_div', 'ks') % average No. of division rounds
plot([0:3], LineFit_all.p1.*[0:3] + LineFit_all.p2 ,'-r') % best fit (average division rate)
xlim([-0.15 3.15]); ylim([0 7.5]);
ylabel('<No. division rounds>'); xlabel('Time (weeks)')

%% ANALYSE SP MODEL GOODNESS-OF-FIT & CELL-CYCLE PERIOD DISTRIBUTION:

%% RUN ABC REJECTION METHOD to deduce progenitor cell-cycle period distributions compatible with the pattern of cell divisions:
if runABC == 1
    % Initialization parameters:
    N = 10; % No. of acceptable posterior estimates (let ~2h calculation for N=100)
    lambda = Lambda_avg; % Range of possible values for the division rate (fixed)
    tlag_range = [0:0.2:2]./7; % Range of possible values for the refractory period (minimum cell-cycle period, in weeks)
    GamShape_range = 2.^[0:0.25:2]; % Range of possible values for the 'Shape' parameter of the gamma-distributed cell-cycle period
    TolThresPreset = [0.0997]; % Tolerance threshold for parameter acceptance/rejection. Preset values are those for ~5% acceptance rate.
    if set_TolThres == 0 % uses predefined optimal tolerance threshold value
        TolThres = TolThresPreset; % Tolerance threshold selection
    else % calculates an adequate tolerance threshold value for 5% acceptance rate
        trialN = 100; % number of trials (random tcc parameter combinations)
        pdist_all = ABCrejection_set_TolThres_acceptRate_BasalCellProlif(rtime,FreqDiv_BL_mean,trialN,lambda,tlag_range,GamShape_range,xtickpos);
        TolThres = quantile(pdist_all,0.05); % sets 5% acceptance rate
    end

    % Run algorithm:
    [OK_lambda,OK_tlag,OK_GamShape] = ABCrejection_tccDist_inference_BasalCellProlif(rtime,FreqDiv_BL_mean,N,TolThres,lambda,tlag_range,GamShape_range,xtickpos);

end

%% SELECT AN OPTIMAL SP-MODEL CELL-CYCLE PERIOD DISTRIBUTION:
if runABC == 1
    % select one of the accepted SP-model cell-cycle period distributions (random choice):
    rndElement = find(mnrnd(1,1/length(OK_tlag)*ones(1,length(OK_tlag)))==1);
    mytlag = OK_tlag(rndElement);
    myGam = OK_GamShape(rndElement);
else % select a pre-specified optimal SP-model cell-cycle period distribution:
    [ParamVal] = SelectModelParamVal('ABC_inference',5); % retrieve optimal tcc parameter values
    Lambda_avg = ParamVal.lambda; mytlag = ParamVal.tlag; myGam = ParamVal.GamShape;
end

%% PLOT OPTIMAL SP-MODEL GAM-tcc FITS ON EXPERIMENTAL H2BGFP-DECONVOLUTED HISTOGRAMS OF No. OF CELL DIVISION ROUNDS:
% Initialization parameters:
M = 10000; % No. of simulated individual basal cells per trial
lambda = Lambda_avg;
tlag = mytlag; % optimal GAM tcc refractory period
GamShape = myGam; % optimal GAM tcc Gamma-shape parameter

% Simulation of SP-model time course in the No. of division rounds:
[NoDiv_BL] = MonteCarloSimulator_SP_BasalCellProlif(rtime,lambda,M,tlag,GamShape);

% actual division values are translated into freqs (FreqDiv_SPmodel) to be compared against experimental freqs (FreqDiv_experim)
myspan = [xtickpos-1 xtickpos(end)];
for buc = 1:length(rtime)
    [counts,centres] = hist(NoDiv_BL(buc,:),myspan);
    FreqDiv_SPmodel{1,buc} = counts(1:end-1) ./ sum(counts(1:end-1));
end

% Plot SP-model fits on H2BGFP-deconvoluted histograms of No. of rounds of cell division:
figure(1)
for buc = 1:length(rtime)
    subplot(3,3,buc)
    hold on
    plot(FreqDiv_SPmodel{1,buc},'-r');
end
