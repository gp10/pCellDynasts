%% SCRIPT USED FOR EVALUATING MODEL GOODNESS-OF-FIT ON SANCHEZ-DANES ET AL (2016) INTERSCALE KRT14 LABELLING PERSISTENCE DATA:
% FITS ON CLONE DENSITY AND FRACTION OF LABELLED BASAL CELLS OVER TIME ARE PLOTTED

%% LOAD EXPERIMENTAL DATA SET ON CLONE DENSITIES:
load ./Datasets/LineageTracing_TailInterscale_K14ClonalDens_SanchezDanes2016_dataset.mat
% ClonalDens: structure array containing variables from the experiment on clone densities (obtained from Sánchez-Danés et al 2016 online publication material)

%% COMPUTATIONAL CLONE SIMULATIONS UNDER SPECIFIC PARAMETER CONDITIONS
% General simulation parameters:
indiv = 10000; % estimated calculation time: ~1min for 10k clones | ~10min for 100k clones
ntime = [0 10.^[-0.8:0.02:1.7]]; %a wide range of time points for which to collect simulation outcome, for best visualization purposes

% The stem cell-committed progenitor (SC-CP) model is simulated under
% specific parameter conditions (those values given by Sánchez-Danés et al
% (2016) and computational clone size distributions retrieved.
ParamSet = SelectModelParamVal('Original',8);
[nx_basal_SCCP,ntime_SCCP] = MonteCarloSimulator_SCCP_BasalCloneDynamics(ntime,ParamSet.labelSfrac,ParamSet.lambdaS,ParamSet.rS,ParamSet.lambdaP,ParamSet.rP,ParamSet.DeltaP,ParamSet.gamma,indiv,ParamSet.tlagS,ParamSet.GamShapeS,ParamSet.tlagP,ParamSet.GamShapeP);

% The single-progenitor (SP) model is simulated under specific parameter
% conditions (MLE values obtained from fitting experimental Krt14 clone size
% distributions in Interscale) and computational clone size distributions
% retrieved.
ParamSet = SelectModelParamVal('MLE',8);
[nx_basal_SP,ntime_SP] = MonteCarloSimulator_SP_BasalCloneDynamics(ntime,ParamSet.lambda,ParamSet.r,ParamSet.gamma,indiv,ParamSet.tlag,ParamSet.GamShape);

%% PLOT FITS ON CLONAL DENSITY OVER TIME, RETRIEVING STATISTICS FROM EXPERIMENTAL DATA:
figure()

% Plot experimental data on clone density:
hold on
all_NClones = [];
all_rtime = [];
for aja = 1:size(ClonalDens.rtime,2)
    % retrieve statistics:
    mean_DensClones(1,aja) = mean(ClonalDens.NClones_I_K14{1,aja});
    std_DensClones(1,aja) = std(ClonalDens.NClones_I_K14{1,aja},0,1);
    sem_DensClones(1,aja) = std_DensClones(1,aja) ./ sqrt(ClonalDens.nmice(aja));
    % plotting:
    plot(repmat(ClonalDens.rtime(aja),ClonalDens.nmice(aja),1),ClonalDens.NClones_I_K14{1,aja},'bo','MarkerFaceColor','b','MarkerSize',3)
    h = errorbar(ClonalDens.rtime(aja),mean_DensClones(1,aja),sem_DensClones(1,aja),'CapSize',0); set(h,'Marker','+','Color','k','LineWidth',0.25);
    %store all values as a vector (for later convenience)
    all_NClones = [all_NClones ClonalDens.NClones_I_K14{1,aja}'];
    all_rtime = [all_rtime repmat(ClonalDens.rtime(aja),1,ClonalDens.nmice(aja))];
end
xlabel('Time (weeks)')
ylabel('Number of clones / mm^2')

% Blind fit, used to retrieve "initial scaling factor" for model fits ("No" will represent the initial population at time 0, which is unknown)
myfun = fittype('No/(1+rlambda*t)','dependent',{'y'},'independent',{'t'},'coefficients',{'No', 'rlambda'});
[myfun_fit,myfun_fit_stat] = fit(all_rtime',all_NClones',myfun,'Startpoint',[80 0.15],'Lower',[0 0],'Upper',[Inf 100]);
N_t0 = myfun_fit.No;

% Plot computational clone density (using SC-CP model):
simClonalDens_SCCP = [];
for aja = 1:size(ntime_SCCP,2)
    simClonalDens_SCCP(1,aja) = size(find(nx_basal_SCCP(:,aja) ~= 0),1);
end
hold on
plot(ntime_SCCP,simClonalDens_SCCP ./ simClonalDens_SCCP(1,1) .* N_t0,'-','Color','g')

% Plot computational clone density (using SP model):
simClonalDens_SP = [];
for aja = 1:size(ntime_SP,2)
    simClonalDens_SP(1,aja) = size(find(nx_basal_SP(:,aja) ~= 0),1);
end
hold on
plot(ntime_SP,simClonalDens_SP ./ simClonalDens_SP(1,1) .* N_t0,'-','Color','r')

xlim([-0.5 25.5]);
ylim([0 70])

%% LOAD EXPERIMENTAL DATA SET ON CLONE SIZES:
% Data from Krt14-driven labelling in interscale  (obtained from Sánchez-Danés et al 2016 online publication material):
load ./Datasets/LineageTracing_TailInterscale_SanchezDanes2016_dataset.mat
% rtime: vector containing experimental timepoints (expressed in weeks)
% nmice: vector containing the number of animals collected per time point
% rx_basal: cell array of format {nmice,timepoints} containing basal clone sizes per individual animal, individual time point
% rx_basal_all: cell array of format {1,timepoints} containing basal clone sizes per individual time point (animals pooled)

%% PLOT FITS ON AVERAGE BASAL CLONE SIZE, RETRIEVING STATISTICS FROM EXPERIMENTAL DATA:
figure()

% Plot experimental data on avg. clone size and retrieve descriptive statistics:
hold on
mean_pooled_CloneSize = [];
sem_pooled_CloneSize = [];
perindiv_mean_CloneSize = {};
mean_indiv_CloneSize = [];
sem_indiv_BCloneSize = [];
for ata = 1:size(rtime,2)
    % average/sem from all clones (pooling from different animals) - most conservative option when considering the experimental error (the one chosen below for error propagation analysis)
    pooled_persis = find(rx_basal_all{1,ata} ~= 0); %only persisting clones (with at least 1 basal cell) are considered
    mean_pooled_CloneSize(ata) = mean(rx_basal_all{1,ata}(pooled_persis,1));
    sem_pooled_CloneSize(ata) = std(rx_basal_all{1,ata}(pooled_persis,1),0) ./ sqrt(size(rx_basal_all{1,ata}(pooled_persis,1),1));
    % average/sem considering data from each animal separately
    for ete = 1:nmice(ata)
        indiv_persis = find(rx_basal{ete,ata} ~= 0); %only persisting clones (with at least 1 basal cell) are considered
        perindiv_mean_CloneSize{1,ata}(ete,1) = mean(rx_basal{ete,ata}(indiv_persis,1));
    end
    mean_indiv_CloneSize(1,ata) = mean(perindiv_mean_CloneSize{1,ata});
    sem_indiv_BCloneSize(1,ata) = std(perindiv_mean_CloneSize{1,ata},0,1) ./ sqrt(nmice(ata));
end
errorbar(rtime,mean_indiv_CloneSize,sem_indiv_BCloneSize,'o','Color','k','MarkerFaceColor','k');
xlabel('Time (weeks)'); ylabel('Avg. basal clone size')

% Plot computational avg. clone size (using SC-CP model):
meanCloneSize_SCCP = [];
for aja = 1:size(ntime_SCCP,2)
    meanCloneSize_SCCP(1,aja) = mean(nx_basal_SCCP(find(nx_basal_SCCP(:,aja) ~= 0),aja));
end
hold on
plot(ntime_SCCP,meanCloneSize_SCCP,'-','Color','g')

% Plot computational avg. clone size (using SP model):
meanCloneSize_SP = [];
for aja = 1:size(ntime_SP,2)
    meanCloneSize_SP(1,aja) = mean(nx_basal_SP(find(nx_basal_SP(:,aja) ~= 0),aja));
end
hold on
plot(ntime_SP,meanCloneSize_SP,'-','Color','r')

xlim([-0.5 25.5]);
ylim([0 8])

%% INFERRED TIME COURSE OF THE %LABELLED CELLS:
% Any trend in the percentage of labelled cells is inferred in the following way:
% %LABELLED CELLS = AvgCloneSize x NoClones (addressing the corresponding error propagation appropriately)
figure()

% Inferred experimental trends:
[NCells_mean,NCells_std] = infer_TotalLabelledCells_trend(mean_pooled_CloneSize,sem_pooled_CloneSize,mean_DensClones,std_DensClones);
errorbar(rtime,NCells_mean./NCells_mean(1),NCells_std./NCells_mean(1),'ko','MarkerFaceColor','k')

% Inferred computational trends (using SP model):
ALL_sampled_NCells_mean_SP = [];
ALL_sampled_NCells_std_SP = [];
for aja = 1:length(rtime); posCollected(aja) = find(ntime' >= rtime(aja)',1); end % retrieve positions in the simulation time point vector matching experimental ones
Nrep = 2000;
for buc = 1:Nrep % repeats Nrep times the subsampling
    
    % Simulate sampling error on average basal clone size:
    NTrackedClones = [];
    for bac = 1:size(rx_basal_all,2) % reproduce the number of clones counted in the experiments
        NTrackedClones = [NTrackedClones size(rx_basal_all{1,bac},1)];
    end
    [sampled_avgCloneSize_SP, sampled_semCloneSize_SP] = simulate_subsampled_avgCloneSize_SanchezDanes2016(nx_basal_SP(:,posCollected),NTrackedClones);
    
    % Simulate sampling error on clone density values considering variable initial labelling induction:
    meanInducedClones_t0 = mean(ClonalDens.NClones_I_K14{1,1});
    stdInducedClones_t0 = std(ClonalDens.NClones_I_K14{1,1},0);
    [sampled_avgNClones_SP, sampled_stdNClones_SP] = simulate_subsampled_InducedCloneDens_SanchezDanes2016(nx_basal_SP(:,posCollected),ClonalDens.nmice,meanInducedClones_t0,stdInducedClones_t0);
    
    % Inference on number of labelled cells (with error propagation):
    [sampled_NCells_mean_SP,sampled_NCells_std_SP] = infer_TotalLabelledCells_trend(sampled_avgCloneSize_SP,sampled_semCloneSize_SP,sampled_avgNClones_SP,sampled_stdNClones_SP);
    
    ALL_sampled_NCells_mean_SP(buc,:) = sampled_NCells_mean_SP;
    ALL_sampled_NCells_std_SP(buc,:) = sampled_NCells_std_SP;
    
end
hold on
% average predicted behavior with the SP
plot(rtime,quantile(ALL_sampled_NCells_mean_SP ./ ALL_sampled_NCells_mean_SP(:,1),0.5,1),'-','Color','r')
% 95% bounds on average predicted behavior with the SP
plot(rtime,quantile(ALL_sampled_NCells_mean_SP ./ ALL_sampled_NCells_mean_SP(:,1),0.025,1),'-','Color','r')
plot(rtime,quantile(ALL_sampled_NCells_mean_SP ./ ALL_sampled_NCells_mean_SP(:,1),0.975,1),'-','Color','r')
% standard margins of uncertainty of the 95% bounds on average predicted behavior with the SP
plot(rtime,quantile(ALL_sampled_NCells_mean_SP ./ ALL_sampled_NCells_mean_SP(:,1),0.025,1) - (mean(ALL_sampled_NCells_std_SP,1) ./ mean(ALL_sampled_NCells_mean_SP,1)),'--','Color','r')
plot(rtime,quantile(ALL_sampled_NCells_mean_SP ./ ALL_sampled_NCells_mean_SP(:,1),0.975,1) + (mean(ALL_sampled_NCells_std_SP,1) ./ mean(ALL_sampled_NCells_mean_SP,1)),'--','Color','r')

ylim([0 4])
ylabel('Labelled cell fraction')
xlabel('Time (weeks)')
xlim([-0.5 25.5])
