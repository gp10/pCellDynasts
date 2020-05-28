%% SCRIPT USED FOR tcc PARAMETER INFERENCE AND FOR SP MODEL GOODNESS-OF-FIT ON H2BGFP DISTRIBUTIONS:
% The average division rate is calculated and an Approximated Bayesian
% Computation (ABC) algorithm can be run to infer the distribution of
% individual progenitor tcc, displaying fits on experimental H2BGFP
% intensity distributions.

%% SELECTION OF EXPERIMENTAL DATA SET AND GENERAL SETTINGS:
selectDataSet = 1; % ( 1=OE_R26rtTA | 2=Paw_R26rtTA | 3=Ear_R26rtTA | 4=Back_R26rtTA | 5=TailScale_R26rtTA | 6=TailInterscale_R26rtTA )
runABC = 0; % ( 1=runs ABC algorithm to compute likeliest cell-cycle period distributions matching H2BGFP dilution profile | 0=skips ABC calculations and loads already preset optimal cell-cycle period distributions referred to in the text )
set_TolThres = 1; % ( 1=Run a test to set an adequate ABC tolerance threshold for a given acceptance rate | 0=Use preset adequate tolerance threshold levels )

%% LOAD EXPERIMENTAL DATA OF H2BGFP INTENSITIES IN BASAL CELL NUCLEI AT DIFFERENT TIMES:
switch selectDataSet
    case 1;        load ./Datasets/H2BGFPdil_OE_R26rtTA_dataset.mat
    case 2;        load ./Datasets/H2BGFPdil_Paw_R26rtTA_dataset.mat
    case 3;        load ./Datasets/H2BGFPdil_Ear_R26rtTA_dataset.mat
    case 4;        load ./Datasets/H2BGFPdil_Back_R26rtTA_dataset.mat
    case 5;        load ./Datasets/H2BGFPdil_TailScale_R26rtTA_dataset.mat
    case 6;        load ./Datasets/H2BGFPdil_TailInterscale_R26rtTA_dataset.mat
end

%% LOG-2 TRANSFORM EXPERIMENTAL DATA, NORMALIZE IT, AND CALCULATE AVERAGE DIVISION RATE:
% Log-2 transformation:
log2_H2BGFP = {};
mean_log2_H2BGFP = {};
for aja = 1:length(rtime)
    mean_log2_H2BGFP_t = [];
    for eje = 1:nmice(aja)
        log2_H2BGFP{eje,aja} = log2(H2BGFP{eje,aja});
        mean_log2_H2BGFP_t = [mean_log2_H2BGFP_t mean(log2_H2BGFP{eje,aja})];
    end
    mean_log2_H2BGFP{1,aja} = mean_log2_H2BGFP_t;
end

% normalize distributions according to average H2BGFP intensity at time0 (done individually per experiment in the case of back-skin data):
nlog2_H2BGFP = {};
mean_nlog2_H2BGFP = {};
nlog2_H2BGFP_all = {};
for aja = 1:length(rtime)
    mean_nlog2_H2BGFP_t = [];
    nlog2_H2BGFP_t = [];
    for eje = 1:nmice(aja)
        if selectDataSet == 4 % back-skin data coming from 3 independent experiments
            if aja <= 3 % for data belonging to experiment A
                nlog2_H2BGFP{eje,aja} = log2_H2BGFP{eje,aja} - mean(mean_log2_H2BGFP{1,1});
            elseif (4 <= aja && aja <= 5) % for data belonging to experiment B
                nlog2_H2BGFP{eje,aja} = log2_H2BGFP{eje,aja} - mean(mean_log2_H2BGFP{1,4});
            else % for data belonging to experiment C
                nlog2_H2BGFP{eje,aja} = log2_H2BGFP{eje,aja} - mean(mean_log2_H2BGFP{1,6});
            end
        else % any other data coming from a single common experiment
            nlog2_H2BGFP{eje,aja} = log2_H2BGFP{eje,aja} - mean(mean_log2_H2BGFP{1,1});
        end
        mean_nlog2_H2BGFP_t = [mean_nlog2_H2BGFP_t mean(nlog2_H2BGFP{eje,aja})];
        nlog2_H2BGFP_t = [nlog2_H2BGFP_t; nlog2_H2BGFP{eje,aja}];
    end
    mean_nlog2_H2BGFP{1,aja} = mean_nlog2_H2BGFP_t;
    nlog2_H2BGFP_all{1,aja} = nlog2_H2BGFP_t;
end

% reshape data format (back-skin data only, so as to compact data from different individual experiments)
if selectDataSet == 4
    [rtime,nmice,H2BGFP,FieldView,log2_H2BGFP,nlog2_H2BGFP,mean_log2_H2BGFP,mean_nlog2_H2BGFP,H2BGFP_all,nlog2_H2BGFP_all] = reshape_BackskinDataFormat(rtime,nmice,H2BGFP,FieldView,log2_H2BGFP,nlog2_H2BGFP,mean_log2_H2BGFP,mean_nlog2_H2BGFP,H2BGFP_all,nlog2_H2BGFP_all);
end

% recenter according to estimated population average at each time point: (only the relative distribution dispersion around the average is accounted for)
rnlog2_H2BGFP = {};
mean_rnlog2_H2BGFP = {};
rnlog2_H2BGFP_all = {};
for aja = 1:length(rtime)
    mean_rnlog2_H2BGFP_t = [];
    rnlog2_H2BGFP_t = [];
    for eje = 1:nmice(aja)
        rnlog2_H2BGFP{eje,aja} = nlog2_H2BGFP{eje,aja} - mean(nlog2_H2BGFP{eje,aja}) + mean(mean_nlog2_H2BGFP{1,aja});
        mean_rnlog2_H2BGFP_t = [mean_rnlog2_H2BGFP_t mean(rnlog2_H2BGFP{eje,aja})];
        rnlog2_H2BGFP_t = [rnlog2_H2BGFP_t; rnlog2_H2BGFP{eje,aja}];
    end
    mean_rnlog2_H2BGFP{1,aja} = mean_rnlog2_H2BGFP_t;
    rnlog2_H2BGFP_all{1,aja} = rnlog2_H2BGFP_t;
end

% Calculate average division rate (from the linear slope of the log2(H2BGFP) decay)
myfun = fittype('p1*x + p2','dependent',{'y'},'independent',{'x'},'coefficients',{'p1', 'p2'});
[LineFit_all,LineFit_all_stat] = fit(repelem(rtime,nmice)',cell2mat(mean_nlog2_H2BGFP)',myfun,'Startpoint',[0 0],'Lower',[-Inf -Inf],'Upper',[Inf Inf]);
lambda_avg = abs(LineFit_all.p1);
disp(sprintf('Avg. Division rate: lambda = %.2f/week',lambda_avg));

% Plot best fit for the average division rate:
figure()
hold on
plot(repelem(rtime,nmice), cell2mat(mean_nlog2_H2BGFP), 'ks') % average H2BGFP intensity
plot([0:3], LineFit_all.p1.*[0:3] + LineFit_all.p2 ,'-r') % best fit (average division rate)
xlim([-0.2 3.2]); ylim([-12 2]); grid on
ylabel('log_2(H2BGFP int.)'); xlabel('Time (weeks)')

%% ANALYSIS OF CELL-CYCLE PERIOD DISTRIBUTION:
if selectDataSet > 4 % tail-skin data sets
    disp('Deconvolution of tcc distribution is not possible for tail skin given large inter-animal variation')
elseif selectDataSet <= 4 % other data sets
    %% RUN ABC REJECTION METHOD to deduce progenitor cell-cycle period distributions compatible with the actual pattern of H2BGFP dilution:
    if runABC == 1
        % Initialization parameters:
        N = 200; % No. of acceptable posterior estimates (let ~2h calculation for N=100)
        lambda = lambda_avg; % Range of possible values for the division rate (fixed)
        tlag_range = [0:0.25:2]./7; % Range of possible values for the refractory period (minimum cell-cycle period, in weeks)
        GamShape_range = 2.^[0:6]; % Range of possible values for the 'Shape' parameter of the gamma-distributed cell-cycle period
        TolThresPreset = [0.159 0.09 0.41 0.36]; % Tolerance threshold for parameter acceptance/rejection. Preset values are those for ~5% acceptance rate.
        if set_TolThres == 0 % uses predefined optimal tolerance threshold value
            TolThres = TolThresPreset(selectDataSet); % Tolerance threshold selection
        else % calculates an adequate tolerance threshold value for 5% acceptance rate
            trialN = 100; % number of trials (random tcc parameter combinations)
            pdist_all = ABCrejection_set_TolThres_acceptRate(rtime,rnlog2_H2BGFP_all,trialN,lambda,tlag_range,GamShape_range);
            TolThres = quantile(pdist_all,0.95); % sets 5% acceptance rate
        end

        % Run algorithm:
        [OK_lambda,OK_tlag,OK_GamShape] = ABCrejection_tccDist_inference(rtime,rnlog2_H2BGFP_all,N,TolThres,lambda,tlag_range,GamShape_range);

        % Store accepted parameter combinations in matrix format:
        OKmatrix = zeros((size(tlag_range,2)-1),size(GamShape_range,2));
        for aja = 1:N
            ycoord = find(GamShape_range == OK_GamShape(aja));
            xcoord = find(OK_tlag(aja) < tlag_range,1);
            OKmatrix(xcoord-1,ycoord) = OKmatrix(xcoord-1,ycoord) + 1;
        end

        % Plot heatmap of most-likely tcc properties (accepted RefractPeriod & GamShape parameter combinations):
        figure()
        imagesc(OKmatrix)
        mycolormap = colormap('gray'); colormap(1-mycolormap)
        set(gca,'XTick',1:7); set(gca,'YTick',1:8);
        set(gca,'XTickLabel',num2cell(GamShape_range)); set(gca,'YTickLabel',{'0.0 - 0.25','0.25 - 0.5','0.5 - 0.75','0.75 - 1.0','1.0 - 1.25','1.25 - 1.5','1.5 - 1.75','1.75 - 2.0'});
        ylabel('Refractory period (days)'); xlabel('Cell-cycle \Gamma-distribution shape, k');
    end

    %% PLOT ACCEPTED (DECONVOLUTED) CELL-CYCLE PERIOD DISTRIBUTIONS:
    figure()
    FigTimeSpan = [0:0.01:12]; % time span the cell-cycle distribution is evaluated for
    hold on
    if runABC == 1
        % plot a representative range of accepted cell-cycle period distributions fitting experimental observations:
        gampdf_tcc = [];
        lagTime_95perc = [];
        for aja = 1:min([50 N])
            gampdf_tcc(aja,:) = gampdf(FigTimeSpan,OK_GamShape(aja),(1/lambda_avg-OK_tlag(aja))./OK_GamShape(aja)); % delayed gamma-distributed tcc
            plot((FigTimeSpan+OK_tlag(aja)).*7,gampdf_tcc(aja,:)./max(gampdf_tcc(aja,:)),'Color',[0.8 0.8 0.8]);
            lagTime_95perc(aja,1) = FigTimeSpan(1,find(cumsum(gampdf_tcc(aja,:),2)>5,1)); % stores lag time (days) for which 95% cells would divide with a longer tcc period (i.e. 5% cutoff in cumulative tcc distribution)
        end
        % select a conservative accepted cell-cycle period distribution:
        [sorted_lagTime_95perc,sorted_lagTime_95perc_I] = sort(lagTime_95perc);
        myloc95perc = sorted_lagTime_95perc_I(find(sorted_lagTime_95perc > quantile(lagTime_95perc,0.05),1)); % selects a distribution within the ~5% of those accepted showing a shortest lag time (widest tcc shape). 
        mytlag = OK_tlag(myloc95perc); myGam = OK_GamShape(myloc95perc);
        mygampdf_tcc = gampdf_tcc(myloc95perc,:); % delayed gamma-distributed tcc
        plot((FigTimeSpan+mytlag).*7,mygampdf_tcc./max(mygampdf_tcc),'Color','r');
    else
        % display the cell-cycle period distribution for pre-specified parameter values:
        [ParamVal] = SelectModelParamVal('ABC_inference',selectDataSet); % retrieve optimal tcc parameter values
        lambda_avg = ParamVal.lambda; mytlag = ParamVal.tlag; myGam = ParamVal.GamShape;
        mygampdf_tcc = gampdf(FigTimeSpan,myGam,(1/lambda_avg-mytlag)./myGam); % delayed gamma-distributed tcc
        plot((FigTimeSpan+mytlag).*7,mygampdf_tcc./max(mygampdf_tcc),'Color','r');
    end
    disp(sprintf('Optimal tcc: lag-time = %.1f days | Gamma-shape = %d',mytlag*7,myGam));

    % show average cell-cycle period time:
    line([1/lambda_avg 1/lambda_avg].*7,[0 1],'Color','b');
    
    hold off
    ylabel('Frequency'); xlabel('Cell-cycle time (days)')
    xlim([0 12]); set(gca,'XTick',[0:2:12]); set(gca,'YTick',[0:0.2:1]);

    %% PLOT OPTIMAL SINGLE PROGENITOR GAM-tcc FITS ON EXPERIMENTAL H2BGFP DILUTION HISTOGRAMS AND COMPARE TO STANDARD EXP-tcc PREDICTION FITS:
    % Initialization parameters:
    M = 10000; % No. of simulated individual basal cells per trial
    tlag2plot = [0, mytlag];% Default EXP vs. optimal GAM tcc refractory period
    GamShape2plot = [1 myGam];% Default EXP vs. optimal GAM tcc Gamma-shape parameter
    mycol = [0.3 0.3 0.3;... % grey color for EXP tcc fit
        1 0 0]; % red color for GAM tcc fit

    figure()
    for sj = 1:length(tlag2plot)

        lambda = lambda_avg;
        tlag = tlag2plot(sj);
        GamShape = GamShape2plot(sj);

        % Simulation of H2BGFP dilution time course under each parameter scenario:
        Idist = MonteCarloSimulator_SP_BasalCellProlif_H2BGFPdil(rtime,rnlog2_H2BGFP_all{1,1},lambda,M,tlag,GamShape);

        % Simulated H2BGFP intensity values are translated into log2 values and these normalized to avg at time0 (as in experimental data)
        clog2_Idist = log2(Idist)-median(log2(Idist(1,:)),2); % simulated
        rlog2_Idist = {}; % experimental
        for buc = 1:length(rtime)
            rlog2_Idist{1,buc} = rnlog2_H2BGFP_all{1,buc}-median(rnlog2_H2BGFP_all{1,1},1);
        end

        % Plot fits on H2BGFP intensity distributions at each time point:
        xrange = [-14:0.25:2];
        for bb = 1:length(rtime) %for each time point
            [simcounts,simcenters] = hist(clog2_Idist(bb,:),xrange); [datacounts,datacenters] = hist(rlog2_Idist{1,bb},xrange);
            subplot(length(rtime),length(rtime),bb)
            hold on
            area(datacenters,datacounts./max(datacounts),'FaceColor','g','EdgeColor','none');
            plot(simcenters,simcounts./max(simcounts),'Color',mycol(sj,:)); 
            hold off
            xlim([-14 2]); ylim([0 1.1])
            for aja = -1:13; line([-aja -aja],[0 1.1],'Color',[0.8 0.8 0.8],'LineStyle','-'); end
            set(gca,'YTick',[0:0.2:1]); set(gca,'XTick',[-10 -5 0])
            box on
            ylabel('Rel. frequency'); xlabel('log_2(H2BGFP int.)')
        end

    end
end
