%% SCRIPT USED FOR EVALUATING MODEL GOODNESS-OF-FIT ON GIROUX ET AL (2017) DATA SET:
% FITS ON CLONE DENSITY AND GFP-CELL DISTRIBUTION ACROSS DIFFERENT LAYERS ARE PLOTTED

%% LOAD EXPERIMENTAL DATA SET:
load ./Datasets/LineageTracing_OE_Giroux2017_dataset.mat
% ClonalDens: structure array containing variables from the experiment on clone densities (digitalization of Fig.S3B from Giroux et al 2017)
% Labelling: structure array containing variables from the experiment on GFP-cell localization (digitalization of Fig.2E from Giroux et al 2017)

%% COMPUTATIONAL CLONE SIMULATIONS UNDER SPECIFIC PARAMETER CONDITIONS
% General simulation parameters:
indiv = 10000; % estimated calculation time: ~1min for 10k clones | ~10min for 100k clones
ntime = [0 10.^[-0.8:0.02:1.7]]; %a wide range of time points for which to collect simulation outcome, for best visualization purposes

% Specify particular model parameter values: (same as MLE obtained from Lrig1 lineage tracing analysis)
ParamSet = SelectModelParamVal('MLE',1);
rho = ParamSet.gamma / (ParamSet.gamma + ParamSet.lambda); % homeostatic condition on the proportion of progenitor basal cells (rho)
m = 1.25; % SB/B cell ratio (experimentally determined from our esophageal data and in agreement with Doupé et al 2012)
mu = rho*ParamSet.lambda/m; % homeostatic condition on shedding rate (/week)

% The single-progenitor (SP) model is simulated under specific parameter
% conditions and computational clone size and cell localization
% distributions retrieved.
[nx_basal,nx_total,nx,ntime] = MonteCarloSimulator_SP_TotalCloneDynamics_Giroux2017(ntime,ParamSet.lambda,ParamSet.r,ParamSet.gamma,mu,indiv,ParamSet.tlag,ParamSet.GamShape);

%% PLOT FITS ON CLONAL DENSITY OVER TIME:
figure()

% Plot experimental data on clone density:
errorbar(ClonalDens.rtime,ClonalDens.avg, abs(ClonalDens.avg-ClonalDens.dnSD),abs(ClonalDens.upSD-ClonalDens.avg) ,'ko','MarkerFaceColor','k');
xlabel('Days post-recombination')
ylabel('#clones/100 basal cells')
box off; xlim([0 200]); ylim([0 4]); set(gca,'XTick',[0:50:200]); set(gca,'YTick',[0:4]);

% Blind fit, used to retrieve "initial scaling factor" for SP fit ("No" will represent the initial population at time 0, which is unknown)
myfun = fittype('No/(1+rlambda*t)','dependent',{'y'},'independent',{'t'},'coefficients',{'No', 'rlambda'});
[myfun_fit,myfun_fit_stat] = fit(ClonalDens.rtime',ClonalDens.avg',myfun,'Startpoint',[3 0.2],'Lower',[1 0],'Upper',[100 Inf]);
N_t0 = myfun_fit.No

% Plot computational clone density (assuming Krt15-driven labelling in the basal layer targets P- and D-cells indistinctively at the beginning):
simClonalDens = [];
for aja = 1:size(ntime,2)
    simClonalDens(1,aja) = size(find(nx_basal(:,aja) ~= 0),1);
end
hold on
plot(ntime*7,simClonalDens ./ simClonalDens(1,1) .* N_t0,'Color','r','LineStyle',':')

% Plot computational clone density (assuming Krt15-driven labelling in the basal layer is restricted to P-cells at the beginning):
simClonalDens = [];
Pcell_derivedClones = find(nx(:,1,1) ~= 0);
for aja = 1:size(ntime,2)
    simClonalDens(1,aja) = size(find(nx_basal(Pcell_derivedClones,aja) ~= 0),1);
end
hold on
plot(ntime*7,simClonalDens ./ simClonalDens(1,1) .* N_t0,'Color','r')

%% PLOT FITS ON GFP-CELL DISTRIBUTION ACROSS LAYERS:
figure()

% Plot experimental data on GFP-cell localization over time:
b = bar([Labelling.B_layer' Labelling.PB_layer' (Labelling.SB_layer+Labelling.SUP_layer)'],'stacked');
set(b(1,1),'FaceColor',[0.2 0.6 1]); set(b,'EdgeColor','none')
set(gca,'XTick',[1:7]); set(gca,'YTick',[0:20:100]);
set(gca,'XTickLabel',{'D1','D3','D5','D7','D14','D28','D56'});
ylabel('% of GFP-positive cells'); xlabel('Days post-recombination')
legend('Basal','Parabasal','Suprabasal','Location','northoutside')

% Plot computational cell-type distributions over time:
simFractionBcells = sum(nx_basal,1) ./ sum(nx_total,1) .* 100;
AnucleatedPercent = 10; %ca. 10% of GFP cells are anucleated SUP cells in the experimental data (these are ignored given that we circumscribe to nucleated suprabasal cells given the nucleated SB/B cell ratio)
simFractionBcells_nucleated = simFractionBcells - AnucleatedPercent;
hold on
for aja = 1:length(Labelling.rtime); posCollected(aja) = find(ntime*7' >= Labelling.rtime(aja)',1); end % retrieve positions in the simulation time point vector matching experimental ones
plot(simFractionBcells_nucleated(posCollected),'-','Color','r')

% Plot SB/B cell ratio corresponding to homeostatic conditions (it has to be rescaled to accout for anucleated SUP cells):
h = line([0 8],[((1/m / (1+1/m)).*100 - 10) ((1/m / (1+1/m)).*100 - 10)]);
set(h,'Color','r','LineStyle','--');