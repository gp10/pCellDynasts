%% SCRIPT USED FOR EVALUATING MODEL GOODNESS-OF-FIT AND FOR PARAMETER INFERENCE:
% THE log-LIKELIHOOD VALUE IS RETRIEVED &/OR FITS ON CLONE SIZE DISTRIBUTIONS ARE PLOTTED

%% SELECTION OF EXPERIMENTAL DATA SET, MODEL FOR FITTING AND OUTPUT SETTINGS:
selectDataSet = 3; % ( 1=OE_Lrig1 | 2=OE_Doupe2012 | 3=Paw_Lim2013 | 4=Ear_Doupe2010 | 5=Back_Murai2018 | 6=Back_Fullgrabe2015 | 7=TailScale_SD2016 | 8=TailInterscale_SD2016 )
selectModelInferMode = 'Original'; % ( 'MLE'=best-fit SP model params. | 'Original'=model params. in original paper | 'Manual'=params. to be set below )
outputLogLike = 0; % ( 0=Does not show the log-Likelihood value, but a detailed display of the fit | 1=shows the log-Likelihood value, limiting the display of the model fit to experimental time points )

%% LOAD EXPERIMENTAL CLONE SIZE DATA AND CONVERT TO CLONE SIZE DISTRIBUTIONS:
switch selectDataSet
    case 1;        load ./Datasets/LineageTracing_OE_Lrig1_dataset.mat
    case 2;        load ./Datasets/LineageTracing_OE_Doupe2012_dataset.mat
    case 3;        load ./Datasets/LineageTracing_Paw_Lim2013_dataset.mat
    case 4;        load ./Datasets/LineageTracing_Ear_Doupe2010_dataset.mat
    case 5;        load ./Datasets/LineageTracing_Back_Murai2018_dataset.mat
    case 6;        load ./Datasets/LineageTracing_Back_Fullgrabe2015_dataset.mat
    case 7;        load ./Datasets/LineageTracing_TailScale_SanchezDanes2016_dataset.mat
    case 8;        load ./Datasets/LineageTracing_TailInterscale_SanchezDanes2016_dataset.mat
end
% rtime: vector containing experimental timepoints (expressed in weeks)
% rx_basal: cell array of format {nmice,timepoints} containing basal clone sizes per individual animal, individual time point
% rx_basal_all: cell array of format {1,timepoints} containing basal clone sizes per individual time point (animals pooled)

% Collect frequencies for each clone size (No. of basal cells):
% we may exclude clones with less than 2 basal cells, since in the experiments they may come from labelling postmitotic, already differentiating cells
[rfreq_all, rfreq_all_rel] = size2freq(rx_basal_all,rtime,1,rx_basal_all,2); % only clones with at least 2 basal cells are considered

% Present frequencies in bins (increasing in powers of 2):
[rfreq_bin_all, rfreq_bin_all_rel, rbin_label] = size2freqbinned(rfreq_all,rx_basal_all,rtime,1);

%% COMPUTATIONAL SIMULATION OF CLONE SIZES UNDER SPECIFIC PARAMETER CONDITIONS
% General simulation parameters:
indiv = 100000; % estimated calculation time: ~1min for 10k clones | ~10min for 100k clones
if outputLogLike == 1
    ntime = rtime; %force the simulation time points to be the same as experimental ones, for MLE calculation
    if (selectDataSet == 3) && strcmpi('Original',selectModelInferMode); ntime = rtime + 1.5/7; end %force the simulation time points to be 1.5d displaced to match experimental ones accounting for the claimed delay in the initial induction in Lim et al (2013)
else
    ntime = [0 10.^[-0.8:0.02:1.7]]; %a wide range of time points for which to collect simulation outcome, for best visualization purposes
end

% Specify particular model parameter values:
if (strcmpi('MLE',selectModelInferMode) || strcmpi('Original',selectModelInferMode))
    % Load given parameter values:
    ParamSet = SelectModelParamVal(selectModelInferMode,selectDataSet);

elseif (strcmpi('Manual',selectModelInferMode))
    % Example with some random parameter values (e.g. to be iterated)
    ParamSet.modelType = 'SP'
    ParamSet.lambda = 2.9; %(/week)
    ParamSet.r = 0.095;
    ParamSet.gamma = 5.3857; %(/week)
    ParamSet.tlag = 0.5/7; % half a day
    ParamSet.GamShape = 8;
end

% Except otherwise stated, the single-progenitor (SP) model is simulated
% under specific parameter conditions and computational clone size
% distributions retrieved. They will constitute the PDF to be tested
% against experimental datasets
if strcmpi('SP',ParamSet.modelType)
    if selectDataSet ~= 6 % start simulations from initial single-cell clones
        [nx_basal,ntime] = MonteCarloSimulator_SP_BasalCloneDynamics(ntime,ParamSet.lambda,ParamSet.r,ParamSet.gamma,indiv,ParamSet.tlag,ParamSet.GamShape);
    else % start simulations from preformed multiple-cell clones (case of analysis of data from Fullgrabe et al 2015)
        time0offset = 3; % experimental time point from which clones are simulated (e.g. 40d post-induction)
        [nx_basal,ntime] = MonteCarloSimulator_SP_BasalCloneDynamics_Fullgrabe2015(ntime,rx_basal_all{1,time0offset},ParamSet.lambda,ParamSet.r,ParamSet.gamma,indiv,ParamSet.tlag,ParamSet.GamShape);
    end
elseif strcmpi('SC-CP',ParamSet.modelType)
    [nx_basal,ntime] = MonteCarloSimulator_SCCP_BasalCloneDynamics(ntime,ParamSet.labelSfrac,ParamSet.lambdaS,ParamSet.rS,ParamSet.lambdaP,ParamSet.rP,ParamSet.DeltaP,ParamSet.gamma,indiv,ParamSet.tlagS,ParamSet.GamShapeS,ParamSet.tlagP,ParamSet.GamShapeP);
end

% Collect frequencies for each clone size (No. of basal cells):
% we may exclude clones with less than 2 basal cells if we did so with the experimental data sets
[nfreq, nfreq_rel] = size2freq(nx_basal,ntime,2,nx_basal,2);

% Present frequencies in bins (increasing in powers of 2):
[nfreq_bin, nfreq_bin_rel, nbin_label] = size2freqbinned(nfreq,nx_basal,ntime,2);

%% GOODNESS-OF-FIT CALCULATED AS log-LIKELIHOOD VALUE:
% A value of log-Likelihood is obtained for the specific set of parameter
% values simulated above.
% This procedure can be repeated for multiple parameter values (if we were
% to loop on them, in a parameter sweep) to estimate the log-Likelihood
% landscape and infer the combination of values yielding the maximum
% likelihood estimate (MLE)

if outputLogLike == 1
    Lbasal_bin_t = [];
    Lbasal_bin = [];
    try
        %                                                       rfreq      myPDF     timepoints
        [Lbasal_bin_t(1,:), Lbasal_bin(1,1)] = logLike_calc(rfreq_bin_all,nfreq_bin_rel,rtime);
    catch
        disp('Could not compute the LogLikelihood!!! Value assigned NaN')
        Lbasal_bin_t(1,:) = NaN(1,size(rtime,2)); Lbasal_bin(1,1) = NaN;
    end

    disp(sprintf('log-Likelihood value: %f',Lbasal_bin))
end

%% PLOT FITS ON BASAL CLONE SIZE DISTRIBUTIONS OVER TIME:
figure()
subplot(2,1,1)
mycol = [1 0 0; 0.93 0.69 0.13; 0.47 0.67 0.19; 0 0 1; 0.75 0 0.75; zeros(10,3)]; % color palette

% Experimental clone size distributions:
avg_rfreq_bin_all_rel = {}; sem_rfreq_bin_all_rel = {};
for aja = 1:length(rtime)
    % mean frequency of clones of certain size:
    avg_rfreq_bin_all_rel{1,aja} = rfreq_bin_all_rel{1,aja}(3:end,1) ./ sum(rfreq_bin_all_rel{1,aja}(3:end,1));
    % SEM of frequency of clones of certain size
    sem_rfreq_bin_all_rel{1,aja} = sqrt(avg_rfreq_bin_all_rel{1,aja} .* (1-avg_rfreq_bin_all_rel{1,aja})) ./ sqrt(sum(rfreq_bin_all{1,aja}(3:end,1)));
    % Plot as errorbars:
    hold on
    for lops = 1:size(avg_rfreq_bin_all_rel{1,aja},1)
        errorbar(rtime(aja),avg_rfreq_bin_all_rel{1,aja}(lops,1),sem_rfreq_bin_all_rel{1,aja}(lops,1),'o','Color',mycol(lops,:),'MarkerFaceColor',mycol(lops,:))
    end
end

% Computational clone size distributions:
avg_nfreq_bin_rel = [];
avg_nfreq_bin_rel = nfreq_bin_rel(3:end,:) ./ sum(nfreq_bin_rel(3:end,:),1);
% Plot frequency time courses:
if selectDataSet ~= 6 % time is recorded in simulation as in the experiments
    for eje = 1:size(avg_nfreq_bin_rel,1)
        if (selectDataSet ~= 3) || strcmpi('MLE',selectModelInferMode) || (outputLogLike==1)
            plot(ntime,avg_nfreq_bin_rel(eje,:),'-','Color',mycol(eje,:))
        else
            plot(ntime+1.5/7,avg_nfreq_bin_rel(eje,:),'-','Color',mycol(eje,:)) % Original fit in Lim et al (2013) considers a 1.5d delay in the initial induction
        end
    end
else % consider the offset in the reference initial time between experiments and simulations
    for eje = 1:size(avg_nfreq_bin_rel,1)
        plot(ntime+rtime(time0offset),avg_nfreq_bin_rel(eje,:),'-','Color',mycol(eje,:))
    end
end

set(gca,'XScale','log'); box on; ylim([0 1]); %xlim([0.3 100]); 
xlabel('Time (weeks)'); ylabel('Clone frequency (%)')
%legend('2','3-4','5-8','9-16','17-32','33-64','>64')

%% PLOT FITS ON AVERAGE BASAL CLONE SIZE OVER TIME:
%figure()
subplot(2,1,2)

% Retrieve and plot experimental avg. clone size:
hold on
mean_indiv_Bclonesizes = {};
mean_Bclonesizes = [];
sem_Bclonesizes = [];
for ata = 1:size(rtime,2)
    for ete = 1:nmice(ata)
        row_indiv_persis = find(rx_basal{ete,ata} ~= 0); %only persisting clones (with at least 1 basal cell) are considered
        mean_indiv_Bclonesizes{1,ata}(ete,1) = mean(rx_basal{ete,ata}(row_indiv_persis,1));
    end
    mean_Bclonesizes(1,ata) = mean(mean_indiv_Bclonesizes{1,ata});
    sem_Bclonesizes(1,ata) = std(mean_indiv_Bclonesizes{1,ata},0,1) ./ sqrt(nmice(ata));
    %plot(repmat(rtime(ata),nmice(ata),1),mean_indiv_Bclonesizes{1,ata},'bo','MarkerFaceColor','b','MarkerSize',3)
end
h = errorbar(rtime,mean_Bclonesizes,sem_Bclonesizes,'o','Color','k','MarkerFaceColor','k');

% Retrieve and plot computational avg. clone size:
mean_All_Bclonesizes = [];
for ata = 1:size(ntime,2)
    row_All_persis = find(nx_basal(:,ata) ~= 0);
    mean_All_Bclonesizes(ata) = mean(nx_basal(row_All_persis,ata));
end
if selectDataSet ~= 6 % time is recorded in simulation as in the experiments
    if (selectDataSet ~= 3) || strcmpi('MLE',selectModelInferMode) || (outputLogLike==1)
        plot(ntime(1,:),mean_All_Bclonesizes,'-','Color','r');
    else
        plot(ntime(1,:)+1.5/7,mean_All_Bclonesizes,'-','Color','r'); % Original fit in Lim et al (2013) considers a 1.5d delay in the initial induction
    end
else % consider the offset in the reference initial time between experiments and simulations
    plot(ntime(1,:)+rtime(time0offset),mean_All_Bclonesizes,'-','Color','r');
end

xlabel('Time (weeks)'); ylabel('Avg. basal clone size')
%legend('Experimental','Computational')
