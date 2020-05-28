%% SCRIPT USED FOR tcc PARAMETER INFERENCE AND FOR SP MODEL GOODNESS-OF-FIT ON DISTRIBUTIONS OF CELL DIVISIONS FROM SADA ET AL (2016):
% The average division rate is calculated and an Approximated Bayesian
% Computation (ABC) algorithm can be run to infer the distribution of
% individual progenitor tcc that best suits the experimental H2BGFP-
% deconvoluted histograms of the No. of division rounds shown for basal
% and suprabasal spinous-layer cells by Sada et al (2016).

%% SELECTION OF GENERAL SETTINGS:
runABC = 0; % ( 1=runs ABC algorithm to compute likeliest SP cell-cycle period distributions and differentiated-cell waiting times for stratification matching cell division profiles across compartments | 0=skips ABC calculations and loads already preset optimal cell-cycle period distributions and waiting times referred to in the text )
set_TolThres = 1; % ( 1=Run a test to set an adequate ABC tolerance threshold for a given acceptance rate | 0=Use preset adequate tolerance threshold levels )
solver_2xSC_CellProlif = 'ODE'; % ('ODE'=No. of division rounds for 2xSC model solved through numerical integration (ODEs) | 'MonteCarlo'=estimates based on Monte-Carlo simulations starting from basal proliferative cells (adequate after the 3-days time point))

%% LOAD H2BGFP-DECONVOLUTED DATA ON FREQUENCIES OF BASAL AND SPINOUS CELL DIVISION NUMBER AT DIFFERENT TIMES FROM Sada et al (2016):
load ./Datasets/H2BGFPdil_Back_Sada2016_dataset.mat

%% PLOT REPORTED HISTOGRAMS OF No. OF CELL DIVISION ROUNDS:
figure(1)
% preset parameters (as in the original paper):
xtickpos{1,1} = 1:5; xtickpos{1,2} = 1:6; xtickpos{1,3} = 1:8;
xlabelType{1,1} = {'0','1','2','3','4'}; xlabelType{1,2} = {'0','1','2','3','4','5'}; xlabelType{1,3} = {'0','1','2','3','4','5','6','7'};
FreqDiv_BL_mean = {}; FreqDiv_SL_mean = {};
rtime = ntime;

% sketch colors (just for legend-display purposes):
ColPalette = 'rbbcr'; StylePalette = {'--','-','--','-','-'};
subplot(3,3,1); hold on; for palb = 1:5; plot([0 0],[0 0],'Color',ColPalette(palb),'LineStyle',StylePalette{1,palb}); end
legend('SP model fit (Sada et al)','2xSC model fit (Sada et al)','2xSC model fit (Sada et al param.)','2xSC model fit (optim param.)','SP model fit (optim param.)')

for aja = 1:length(ntime) %different time points
    
    % BL
    subplot(3,3,3+aja)
    hold on
    FreqDiv_BL_all = [];
    for aje = 1:nmice_BL(aja)
        FreqDiv_BL_all(xtickpos{1,aja},aje) = FreqDiv_BL{aje,aja};
    end
    FreqDiv_BL_mean{1,aja} = mean(FreqDiv_BL_all,2);
    bar(FreqDiv_BL_mean{1,aja},'FaceColor','g'); ylim([0 0.9])
    for aje = 1:nmice_BL(aja)
        plot(xtickpos{1,aja},FreqDiv_BL{aje,aja},'ko','MarkerFaceColor','k','MarkerSize',3)
    end
    xlim([xtickpos{1,aja}(1)-0.7 xtickpos{1,aja}(end)+0.7]); set(gca,'XTick',xtickpos{1,aja}); set(gca,'XTickLabel',xlabelType{1,aja});
    xlabel('No. cell divisions'); ylabel('Frequency')
    
    % SL
    subplot(3,3,aja)
    hold on
    FreqDiv_SL_all = [];
    for aje = 1:nmice_SL(aja)
        FreqDiv_SL_all(xtickpos{1,aja},aje) = FreqDiv_SL{aje,aja};
    end
    FreqDiv_SL_mean{1,aja} = mean(FreqDiv_SL_all,2);
    bar(FreqDiv_SL_mean{1,aja},'FaceColor','g'); ylim([0 0.9])
    for aje = 1:nmice_SL(aja)
        plot(xtickpos{1,aja},FreqDiv_SL{aje,aja},'ko','MarkerFaceColor','k','MarkerSize',3)
    end
    xlim([xtickpos{1,aja}(1)-0.7 xtickpos{1,aja}(end)+0.7]); set(gca,'XTick',xtickpos{1,aja}); set(gca,'XTickLabel',xlabelType{1,aja});
    xlabel('No. cell divisions'); ylabel('Frequency')
    
end

%% PLOT SADA ET AL (2016) SP MODEL AND 2xSC MODEL FITS AS DISPLAYED IN THE ORIGINAL PUBLICATION:
for aja = 1:length(ntime) % different time points
    
    % BL
    subplot(3,3,3+aja)
    hold on
    plot(FreqDiv_BL_SPmodel{1,aja},'r','LineStyle','--'); % original SP
    plot(FreqDiv_BL_2xSCmodel{1,aja},'b'); % original 2xSC
    
    % SL
    subplot(3,3,aja)
    hold on
    plot(FreqDiv_SL_SPmodel{1,aja},'r','LineStyle','--'); % original SP
    plot(FreqDiv_SL_2xSCmodel{1,aja},'b'); % original 2xSC

end

%% PLOT ACTUAL 2xSC MODEL FIT CORRESPONDING TO THE PARAMETERS GIVEN IN SADA ET AL (2016):

% load 2xSC model parameters (semi-coupled version):
[ParamVal] = SelectModelParamVal('Original/Ad-hoc',5);

if strcmpi('MonteCarlo',solver_2xSC_CellProlif) % Solver based on Monte-Carlo simulations: (approximate; adequate for long-term)
    
    % run simulations of basal cell proliferation according to model parameters:
    indiv = 100000; % number of basal cells to simulate (allow aprox. ~2min calculation per model run)
    [nx_basal,nx_total,ntime,X,Tini,Tend,AllNoDiv] = MonteCarloSimulator_2xSC_TotalCloneDynamics_CellProlif(rtime,ParamVal.labelS1frac,ParamVal.lambdaS1,ParamVal.lambdaS2,ParamVal.uS1,ParamVal.gammaS1,ParamVal.gammaS2,ParamVal.mu,indiv,ParamVal.tlagDivS1,ParamVal.GamShapeDivS1,ParamVal.tlagStrS1,ParamVal.GamShapeStrS1,ParamVal.tlagDivS2,ParamVal.GamShapeDivS2,ParamVal.tlagStrS2,ParamVal.GamShapeStrS2,ParamVal.tlagShed,ParamVal.GamShapeShed);

    % retrieve No. of division rounds in basal and suprabasal cells:
    NoDiv_BL = {}; NoDiv_SL = {};
    for aja = 1:length(ntime)
        NoDiv_BL_mytime = []; NoDiv_SL_mytime = [];
        for lep = 1:size(X,2)
            try
                myTargetCells_BL = find( (Tini{1,lep} <= rtime(1,aja)) & (rtime(1,aja) < Tend{1,lep}) & X{1,lep}(:,3)==0 ); % basal cells
                myTargetCells_SL = find( (Tini{1,lep} <= rtime(1,aja)) & (rtime(1,aja) < Tend{1,lep}) & X{1,lep}(:,3)==1 ); % suprabasal (spinous) cells
                mySubsetDiv_BL = AllNoDiv{1,lep}(myTargetCells_BL,1);
                mySubsetDiv_SL = AllNoDiv{1,lep}(myTargetCells_SL,1);
                NoDiv_BL_mytime = [NoDiv_BL_mytime; mySubsetDiv_BL];
                NoDiv_SL_mytime = [NoDiv_SL_mytime; mySubsetDiv_SL];
            end
        end
        NoDiv_BL{1,aja} = NoDiv_BL_mytime;
        NoDiv_SL{1,aja} = NoDiv_SL_mytime;
    end

    % plot frequency in the No. of division rounds:
    for buc = 1:length(rtime) % different time points

        % BL
        subplot(3,3,3+buc)
        hold on
        % actual division values are translated into freqs (FreqDiv_2xSCmodel) to be compared against experimental freqs (FreqDiv_experim)
        myspan = [xtickpos{1,buc}-1 xtickpos{1,buc}(end)];
        [counts,centres] = hist(NoDiv_BL{1,buc},myspan);
        FreqDiv_BL_2xSCmodel_v2{1,buc} = counts(1:end-1) ./ sum(counts(1:end-1));
        plot(FreqDiv_BL_2xSCmodel_v2{1,buc},'b','LineStyle','--');

        % SL
        subplot(3,3,buc)
        hold on
        % actual division values are translated into freqs (FreqDiv_2xSCmodel) to be compared against experimental freqs (FreqDiv_experim)
        myspan = [xtickpos{1,buc}-1 xtickpos{1,buc}(end)];
        [counts,centres] = hist(NoDiv_SL{1,buc},myspan);
        FreqDiv_SL_2xSCmodel_v2{1,buc} = counts(1:end-1) ./ sum(counts(1:end-1));
        plot(FreqDiv_SL_2xSCmodel_v2{1,buc},'b','LineStyle','--');

    end
    
elseif strcmpi('ODE',solver_2xSC_CellProlif) % Solver based on numerical integration of ODEs: (exact)
    
    % set initial population fractions (No. of divisions) in homeostasis:
    x0_div = zeros(45,1);
    x0_div(1) = 1; %densS1a_div0(t0) = 1 (arbitrary)
    x0_div(10) = x0_div(1); %densS1b_div0(t0) (ratio given at homeostasis)
    x0_div(19) = (1-ParamVal.labelS1frac)/ParamVal.labelS1frac * (x0_div(1)+x0_div(10)); %densS2_div0(t0) (ratio given at homeostasis)
    x0_div(28) = (2*ParamVal.lambdaS1*x0_div(10) + 0.5*ParamVal.gammaS2*x0_div(19)) / (2*ParamVal.mu); %densSLa_div0(t0) (ratio given at homeostasis)
    x0_div(37) = x0_div(28); %densSLb_div0(t0) (ratio given at homeostasis)

    % numerical integration (ODEs) to get No. of division rounds in basal and suprabasal cells:
    myTimeSpan = [0:0.01:6]; %weeks
    ode=@(t,x) ODEs_2xSC_CellProlif_GAMtcc(t,x,ParamVal.lambdaS1,ParamVal.lambdaS2,ParamVal.gammaS2,ParamVal.mu);
    [t,NoDiv]=ode45(ode,myTimeSpan,x0_div);
    
    % plot frequency in the No. of division rounds:
    for buc = 1:length(rtime) % different time points
        
        loc2plot = find(myTimeSpan>rtime(buc),1);

        % BL ( S1a basal + S1b basal + S2 basal )
        subplot(3,3,3+buc)
        hold on
        plot( (NoDiv(loc2plot,xtickpos{1,buc})+NoDiv(loc2plot,xtickpos{1,buc}+9)+NoDiv(loc2plot,xtickpos{1,buc}+18)) ./ sum(NoDiv(loc2plot,xtickpos{1,buc})+NoDiv(loc2plot,xtickpos{1,buc}+9)+NoDiv(loc2plot,xtickpos{1,buc}+18)), 'b', 'LineStyle','--')

        % SL ( SLa spinous + SLb spinous )
        subplot(3,3,buc)
        hold on
        plot( (NoDiv(loc2plot,xtickpos{1,buc}+27)+NoDiv(loc2plot,xtickpos{1,buc}+36)) ./ sum(NoDiv(loc2plot,xtickpos{1,buc}+27)+NoDiv(loc2plot,xtickpos{1,buc}+36)), 'b', 'LineStyle','--')

    end

end

%% PLOT 2xSC MODEL FIT CORRESPONDING TO THE PARAMETERS GIVEN IN SADA ET AL (2016) BUT WITH lambdaS2 TWICE AS LARGE:
% (these parameter conditions should reproduce the 2xSC model fits shown by Sada et al in their original publication)

% load 2xSC model parameters (semi-coupled version):
[ParamVal] = SelectModelParamVal('Original/Ad-hoc',6);

if strcmpi('MonteCarlo',solver_2xSC_CellProlif) % Solver based on Monte-Carlo simulations: (approximate; adequate for long-term)

    % run simulations of basal cell proliferation according to that model:
    indiv = 100000; % number of basal cells to simulate (allow aprox. ~2min calculation per model run)
    [nx_basal,nx_total,ntime,X,Tini,Tend,AllNoDiv] = MonteCarloSimulator_2xSC_TotalCloneDynamics_CellProlif(rtime,ParamVal.labelS1frac,ParamVal.lambdaS1,ParamVal.lambdaS2,ParamVal.uS1,ParamVal.gammaS1,ParamVal.gammaS2,ParamVal.mu,indiv,ParamVal.tlagDivS1,ParamVal.GamShapeDivS1,ParamVal.tlagStrS1,ParamVal.GamShapeStrS1,ParamVal.tlagDivS2,ParamVal.GamShapeDivS2,ParamVal.tlagStrS2,ParamVal.GamShapeStrS2,ParamVal.tlagShed,ParamVal.GamShapeShed);

    % retrieve No. of division rounds in basal and suprabasal cells:
    NoDiv_BL = {}; NoDiv_SL = {};
    for aja = 1:length(ntime)
        NoDiv_BL_mytime = []; NoDiv_SL_mytime = [];
        for lep = 1:size(X,2)
            try
                myTargetCells_BL = find( (Tini{1,lep} <= rtime(1,aja)) & (rtime(1,aja) < Tend{1,lep}) & X{1,lep}(:,3)==0 ); % basal cells
                myTargetCells_SL = find( (Tini{1,lep} <= rtime(1,aja)) & (rtime(1,aja) < Tend{1,lep}) & X{1,lep}(:,3)==1 ); % suprabasal (spinous) cells
                mySubsetDiv_BL = AllNoDiv{1,lep}(myTargetCells_BL,1);
                mySubsetDiv_SL = AllNoDiv{1,lep}(myTargetCells_SL,1);
                NoDiv_BL_mytime = [NoDiv_BL_mytime; mySubsetDiv_BL];
                NoDiv_SL_mytime = [NoDiv_SL_mytime; mySubsetDiv_SL];
            end
        end
        NoDiv_BL{1,aja} = NoDiv_BL_mytime;
        NoDiv_SL{1,aja} = NoDiv_SL_mytime;
    end

    % plot frequency in the No. of division rounds:
    for buc = 1:length(rtime) % different time points

        % BL
        subplot(3,3,3+buc)
        hold on
        % actual division values are translated into freqs (FreqDiv_2xSCmodel) to be compared against experimental freqs (FreqDiv_experim)
        myspan = [xtickpos{1,buc}-1 xtickpos{1,buc}(end)];
        [counts,centres] = hist(NoDiv_BL{1,buc},myspan);
        FreqDiv_BL_2xSCmodel_v3{1,buc} = counts(1:end-1) ./ sum(counts(1:end-1));
        plot(FreqDiv_BL_2xSCmodel_v3{1,buc},'c');

        % SL
        subplot(3,3,buc)
        hold on
        % actual division values are translated into freqs (FreqDiv_2xSCmodel) to be compared against experimental freqs (FreqDiv_experim)
        myspan = [xtickpos{1,buc}-1 xtickpos{1,buc}(end)];
        [counts,centres] = hist(NoDiv_SL{1,buc},myspan);
        FreqDiv_SL_2xSCmodel_v3{1,buc} = counts(1:end-1) ./ sum(counts(1:end-1));
        plot(FreqDiv_SL_2xSCmodel_v3{1,buc},'c');

    end
    
elseif strcmpi('ODE',solver_2xSC_CellProlif) % Solver based on numerical integration of ODEs: (exact)

    % set initial population fractions (No. of divisions) in homeostasis:
    x0_div = zeros(45,1);
    x0_div(1) = 1; %densS1a_div0(t0) = 1 (arbitrary)
    x0_div(10) = x0_div(1); %densS1b_div0(t0) (ratio given at homeostasis)
    x0_div(19) = (1-ParamVal.labelS1frac)/ParamVal.labelS1frac * (x0_div(1)+x0_div(10)); %densS2_div0(t0) (ratio given at homeostasis)
    x0_div(28) = (2*ParamVal.lambdaS1*x0_div(10) + 0.5*ParamVal.gammaS2*x0_div(19)) / (2*ParamVal.mu); %densSLa_div0(t0) (ratio given at homeostasis)
    x0_div(37) = x0_div(28); %densSLb_div0(t0) (ratio given at homeostasis)

    % numerical integration (ODEs) to get No. of division rounds in basal and suprabasal cells:
    myTimeSpan = [0:0.01:6]; %weeks
    ode=@(t,x) ODEs_2xSC_CellProlif_GAMtcc(t,x,ParamVal.lambdaS1,ParamVal.lambdaS2,ParamVal.gammaS2,ParamVal.mu);
    [t,NoDiv]=ode45(ode,myTimeSpan,x0_div);
    
    % plot frequency in the No. of division rounds:
    for buc = 1:length(rtime) % different time points
        
        loc2plot = find(myTimeSpan>rtime(buc),1);

        % BL ( S1a basal + S1b basal + S2 basal )
        subplot(3,3,3+buc)
        hold on
        plot( (NoDiv(loc2plot,xtickpos{1,buc})+NoDiv(loc2plot,xtickpos{1,buc}+9)+NoDiv(loc2plot,xtickpos{1,buc}+18)) ./ sum(NoDiv(loc2plot,xtickpos{1,buc})+NoDiv(loc2plot,xtickpos{1,buc}+9)+NoDiv(loc2plot,xtickpos{1,buc}+18)), 'c')

        % SL ( SLa spinous + SLb spinous )
        subplot(3,3,buc)
        hold on
        plot( (NoDiv(loc2plot,xtickpos{1,buc}+27)+NoDiv(loc2plot,xtickpos{1,buc}+36)) ./ sum(NoDiv(loc2plot,xtickpos{1,buc}+27)+NoDiv(loc2plot,xtickpos{1,buc}+36)), 'c')

    end

end

%% ESTIMATE AVERAGE DIVISION RATE (exclusively from BL data)
% Retrieve average number of divisions over time:
avg_No_div = [];
for aja = 1:length(ntime)
    avg_No_div(1,aja) = sum( FreqDiv_BL_mean{1,aja}.* (xtickpos{1,aja}-1)' ) ./ sum(FreqDiv_BL_mean{1,aja});
end
% Linear fit on the average number of divisions over time (fit restricted to t=3 & t=7, as at 21d the histograms are cut up to 7 divisions)
myfun = fittype('p1*x + p2','dependent',{'y'},'independent',{'x'},'coefficients',{'p1', 'p2'});
[LineFit_all,LineFit_all_stat] = fit(ntime(1:2)',avg_No_div(1:2)',myfun,'Startpoint',[0 0],'Lower',[-Inf 0],'Upper',[Inf 0]);
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

%% RUN ABC REJECTION METHOD to deduce progenitor cell-cycle period distributions compatible with the pattern of cell divisions in basal layer:
if runABC == 1
    % Initialization parameters:
    N = 10; % No. of acceptable posterior estimates (let ~2h calculation for N=100)
    lambda = Lambda_avg; % Range of possible values for the division rate (fixed)
    tlag_range = [0:0.2:2]./7; % Range of possible values for the refractory period (minimum cell-cycle period, in weeks)
    GamShape_range = 2.^[0:0.25:2]; % Range of possible values for the 'Shape' parameter of the gamma-distributed cell-cycle period
    TolThresPreset = [0.0539]; % Tolerance threshold for parameter acceptance/rejection. Preset values are those for ~5% acceptance rate.
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
    [ParamVal] = SelectModelParamVal('ABC_inference',6); % retrieve optimal tcc parameter values
    Lambda_avg = ParamVal.lambda; mytlag = ParamVal.tlag; myGam = ParamVal.GamShape;
end

%% PLOT OPTIMAL SP-MODEL GAM-tcc FITS ON EXPERIMENTAL H2BGFP-DECONVOLUTED HISTOGRAMS OF No. OF CELL DIVISION ROUNDS in BASAL LAYER:
% Initialization parameters:
M = 10000; % No. of simulated individual basal cells per trial
lambda = Lambda_avg;
tlag = mytlag; % optimal GAM tcc refractory period
GamShape = myGam; % optimal GAM tcc Gamma-shape parameter

% Simulation of SP-model time course in the No. of division rounds:
[NoDiv_BL] = MonteCarloSimulator_SP_BasalCellProlif(rtime,lambda,M,tlag,GamShape);

% actual division values are translated into freqs (FreqDiv_SPmodel) to be compared against experimental freqs (FreqDiv_experim)
for buc = 1:length(rtime)
    myspan = [xtickpos{1,buc}-1 xtickpos{1,buc}(end)];
    [counts,centres] = hist(NoDiv_BL(buc,:),myspan);
    FreqDiv_SPmodel{1,buc} = counts(1:end-1) ./ sum(counts(1:end-1));
end

% Plot SP-model fits on H2BGFP-deconvoluted histograms of No. of rounds of cell division:
figure(1)
for buc = 1:length(rtime)
    subplot(3,3,3+buc)
    hold on
    plot(FreqDiv_SPmodel{1,buc},'-r');
end

%% RUN ABC REJECTION METHOD to deduce stratification/shedding waiting time period distributions compatible with the pattern of cell divisions in spinous (suprabasal) layer:
if runABC == 1
    % Initialization parameters:
    N = 10; % No. of acceptable posterior estimates (let ~2h calculation for N=10)
    lambda = Lambda_avg; % optimal average division rate, as obtained above (fixed)
    tlagDiv = mytlag; % optimal GAM tcc refractory period, as obtained above (fixed)
    GamShapeDiv = myGam; % optimal GAM tcc Gamma-shape parameter, as obtained above (fixed)
    mu_range = [2.^[0:1:5]].*Lambda_avg; % Range of possible values for the average shedding rate (/week)
    tlagShed_range = [0:0.3:0.6]; % Range of possible values for the waiting time (minimum period, in weeks) before stratification/shedding, relative to avg. stratification/shedding time (1/mu)
    GamShapeShed_range = 2.^[0:0.5:2]; % Range of possible values for the 'Shape' parameter of the gamma-distributed stratification/shedding waiting time
    TolThresPreset = [0.0557]; % Tolerance threshold for parameter acceptance/rejection. Preset values are those for ~5% acceptance rate.
    if set_TolThres == 0 % uses predefined optimal tolerance threshold value
        TolThres = TolThresPreset; % Tolerance threshold selection
    else % calculates an adequate tolerance threshold value for 5% acceptance rate
        trialN = 100; % number of trials (random waiting-time parameter combinations)
        pdist_all = ABCrejection_set_TolThres_acceptRate_SuprabasalCellDiv(rtime,FreqDiv_SL_mean,trialN,lambda,tlagDiv,GamShapeDiv,mu_range,tlagShed_range,GamShapeShed_range,xtickpos);
        TolThres = quantile(pdist_all,0.05); % sets 5% acceptance rate
    end

    % Run algorithm:
    [OK_mu,OK_tlagShed,OK_GamShapeShed] = ABCrejection_tccDist_inference_SuprabasalCellDiv(rtime,FreqDiv_SL_mean,N,TolThres,lambda,tlagDiv,GamShapeDiv,mu_range,tlagShed_range,GamShapeShed_range,xtickpos);

end

%% SELECT AN OPTIMAL SP-MODEL STRATIFICATION/SHEDDING TIME PERIOD DISTRIBUTION:
if runABC == 1
    % select one of the accepted SP-model waiting time period distributions (random choice):
    rndElement = find(mnrnd(1,1/length(OK_mu)*ones(1,length(OK_mu)))==1);
    mymu = OK_mu(rndElement);
    mytlagShed = OK_tlagShed(rndElement);
    myGamShapeShed = OK_GamShapeShed(rndElement);
else % select a pre-specified optimal SP-model waiting time period distribution:
    [ParamVal] = SelectModelParamVal('ABC_inference',6); % retrieve optimal waiting time parameter values
    mymu = ParamVal.mu; mytlagShed = ParamVal.tlagShed; myGamShapeShed = ParamVal.GamShapeShed;
end

%% PLOT OPTIMAL SP-MODEL GAM-SHEDDING-WAITING-TIME FITS ON EXPERIMENTAL H2BGFP-DECONVOLUTED HISTOGRAMS OF No. OF CELL DIVISION ROUNDS in SPINOUS LAYER:
% Initialization parameters:
M = 1000; % No. of simulated individual spinous cells per trial
lambda = Lambda_avg;
tlagDiv = mytlag; % optimal GAM tcc refractory period
GamShapeDiv = myGam; % optimal GAM tcc Gamma-shape parameter
mu = mymu; % optimal average shedding rate
tlagShed = mytlagShed; % optimal refractory period before shedding
GamShapeShed = myGamShapeShed; % optimal Gamma-shape parameter of shedding waiting time 
rho = mu / (lambda + mu); % homeostatic requirement

% Simulation of SP-model time course in the No. of division rounds in spinous (suprabasal) compartment:
[NoDiv_SL] = MonteCarloSimulator_SP_SuprabasalCellDiv(rtime,M,rho,lambda,tlagDiv,GamShapeDiv,mu,tlagShed,GamShapeShed);

% actual division values are translated into freqs (FreqDiv_SPmodel) to be compared against experimental freqs (FreqDiv_experim)
for buc = 1:length(rtime)
    myspan = [xtickpos{1,buc}-1 xtickpos{1,buc}(end)];
    [counts,centres] = hist(NoDiv_SL(buc,:),myspan);
    FreqDiv_SPmodel{1,buc} = counts(1:end-1) ./ sum(counts(1:end-1));
end

% Plot SP-model fits on H2BGFP-deconvoluted histograms of No. of rounds of cell division in spinous layer:
figure(1)
for buc = 1:length(rtime)
    subplot(3,3,buc)
    hold on
    plot(FreqDiv_SPmodel{1,buc},'-r');
end
