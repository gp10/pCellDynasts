%% SCRIPT USED FOR EXPLORING H2BGFP DILUTION PROFILES UNDER DIFFERENT THEORETICAL MODEL PARADIGMS:
% Cell-proliferation dynamics are simulated under different theoretical
% scenarios, allowing to explore the pattern of H2BGFP dilution distributions
% expected under a SP model vs. SC-CP or 2xSC model paradigms.

%% SELECTION OF GENERAL SETTINGS:
setModel2plot = [1,2,3,4]; % ( 1=SP_model | 2=SC-CP_Mascre2012 | 3=SC-CP_SD2016 | 4=2xSC_Sada2016 )
selectTime2plot = 3; %(weeks) (time of chase at which to plot distribution of individual H2BGFP intensities)
set_rnlog2_H2BGFP_t0 = 0; % ( 1=customized prior distribution of log2 H2BGFP intensities at time=0 passed as a column vector | 0=normalized keratinocyte log2(H2BGFP) intensity values in OE at time0 are taken by default )
indiv = 100000; % number of basal cells to simulate (allow aprox. ~2min calculation per model run)
rtime = [0.2857    0.4286    0.7143    1.0000    1.4286    2.0000    3.0000    4.0000]; %(weeks) (time points for clone size collection (arbitrary here, but they should span beyond the time of chase of interest)

%% RUN H2BGFP PREDICITONS WITH EACH THERETICAL MODEL:
Idist = {};
Idist_SC = {};
Idist_PD = {};
for aja = setModel2plot

    %% REFERENCE H2BGFP DISTRIBUTION USED AS PRIOR FOR TIME=0:
    if set_rnlog2_H2BGFP_t0 == 0
        load('./Datasets/rnlog2_H2BGFP_t0_OE.mat','rnlog2_H2BGFP_t0'); % intensities taken from experimental data from OE at time0
    else
        rnlog2_H2BGFP_t0 = randn(10000,1)./4; % customized distribution of log2 H2BGFP intensity values at time0
    end
    
    %% LOAD AND RUN CELL PROLIFERATION DYNAMICS ACCORDING TO EACH MODEL PARADIGM:
    % load specific model:
    [ParamVal] = SelectModelParamVal('Original/Ad-hoc',aja);
    % run simulations of that model:
    if strcmpi('SP',ParamVal.modelType) % run SP model
        [nx_basal,ntime,X,Tini,Tend,AllNoDiv] = MonteCarloSimulator_SP_BasalCloneDynamics_CellProlif(rtime,ParamVal.lambda,ParamVal.r,ParamVal.gamma,indiv,ParamVal.tlag,ParamVal.GamShape);
    elseif strcmpi('SC-CP',ParamVal.modelType) % run SC-CP model from Mascre2012 (aja=2) or Sanchez-Danes2016 (aja=3)
        [nx_basal,ntime,X,Tini,Tend,AllNoDiv] = MonteCarloSimulator_SCCP_BasalCloneDynamics_CellProlif(rtime,ParamVal.labelSfrac,ParamVal.lambdaS,ParamVal.rS,ParamVal.lambdaP,ParamVal.rP,ParamVal.DeltaP,ParamVal.gamma,indiv,ParamVal.tlagS,ParamVal.GamShapeS,ParamVal.tlagP,ParamVal.GamShapeP);
    elseif strcmpi('2xSC',ParamVal.modelType) % run 2xSC model from Sada2016
        [nx_basal,nx_total,ntime,X,Tini,Tend,AllNoDiv] = MonteCarloSimulator_2xSC_TotalCloneDynamics_CellProlif(rtime,ParamVal.labelS1frac,ParamVal.lambdaS1,ParamVal.lambdaS2,ParamVal.uS1,ParamVal.gammaS1,ParamVal.gammaS2,ParamVal.mu,indiv,ParamVal.tlagDivS1,ParamVal.GamShapeDivS1,ParamVal.tlagStrS1,ParamVal.GamShapeStrS1,ParamVal.tlagDivS2,ParamVal.GamShapeDivS2,ParamVal.tlagStrS2,ParamVal.GamShapeStrS2,ParamVal.tlagShed,ParamVal.GamShapeShed);
    end

    %% RETRIEVE No DIVISION ROUNDS EXPERIENCED BY EACH BASAL CELL DURING CHASE AND SIMULATE ITS CORRESPONDING H2BGFP DECAY PROCESS:
    % parameter definition:
    ALL_It = [];
    ALL_It_SC = [];
    ALL_It_PD = [];
    ALL_type = [];

    for aje = 1:size(X,2) % for each clone
        try
            % retrieve all living basal cell locators
            if (strcmpi('2xSC',ParamVal.modelType) ~= 1)
                myTargetCells = find( (Tini{1,aje} <= selectTime2plot) & (selectTime2plot < Tend{1,aje}) );
            else
                myTargetCells = find( (Tini{1,aje} <= selectTime2plot) & (selectTime2plot < Tend{1,aje}) & X{1,aje}(:,3)==0 );
            end
            
            % retrieve no. of division rounds experienced by each living basal cell
            mySubsetDiv = AllNoDiv{1,aje}(myTargetCells,1);
            
            % simulate H2BGFP dilution process
            for iji = 1:size(mySubsetDiv,1)
                It = rnlog2_H2BGFP_t0(randperm(size(rnlog2_H2BGFP_t0,1),1),1); % pick a random initial log2 H2BGFP intensity
                It = 2.^It; % log2 values are back-translated into actual values
                counter = 0;
                while mySubsetDiv(iji)>counter
                    It = It/2 + (rand-0.5)*sqrt(12)*0.06*It; % two-fold dilution upon every round of cell division
                    counter = counter + 1;
                end
                
                % store final H2BGFP intensity
                ALL_It = [ALL_It It];
                if ((strcmpi('SC-CP',ParamVal.modelType)) && (X{1,aje}(myTargetCells(iji),1)==1)) || ((strcmpi('2xSC',ParamVal.modelType)) && (X{1,aje}(myTargetCells(iji),2)==1))
                    ALL_It_SC = [ALL_It_SC It]; %it's a slow-cycling stem cell (SC)
                    ALL_type = [ALL_type 1];
                else
                    ALL_It_PD = [ALL_It_PD It]; %it's either a quickly-dividing SC or progenitor cell (P) or a differentiated cell (D)
                    ALL_type = [ALL_type 2];
                end
                
            end
        end
    end
    Idist{1,aja} = ALL_It;
    Idist_SC{1,aja} = ALL_It_SC;
    Idist_PD{1,aja} = ALL_It_PD;

    %% RESULTING H2BGFP DISTRIBUTIONS ARE log2 TRANSFORMED & NORMALIZED TO AVERAGE AT TIME=0 (we are just interested in their shape/dispersion)
    cdist = log2(Idist{1,aja})-median(rnlog2_H2BGFP_t0);
    cdist_SC = log2(Idist_SC{1,aja})-median(rnlog2_H2BGFP_t0);
    cdist_PD = log2(Idist_PD{1,aja})-median(rnlog2_H2BGFP_t0);

    %% PLOT RESULTING NORMALIZED H2BGFP DISTRIBUTIONS UNDER THE DIFFERENT MODEL HYPOTHESES:
    xrange = [-14:0.25:2];
    [simcounts,simcenters] = hist(cdist,xrange);
    [simcounts_SC,simcenters_SC] = hist(cdist_SC,xrange);
    [simcounts_PD,simcenters_PD] = hist(cdist_PD,xrange);

    figure(1)
    subplot(length(setModel2plot),length(setModel2plot),aja)
    hold on
    plot(simcenters,simcounts./max(simcounts),'k');
    plot(simcenters_SC,simcounts_SC./max(simcounts),'r'); % any slow-cycling subpopulation
    plot(simcenters_PD,simcounts_PD./max(simcounts),'b'); % any quickly-dividing subpopulation
    hold off
    title(ParamVal.modelType)
    xlim([-14 2]); ylim([0 1.1])
    for sjs = -1:13; line([-sjs -sjs],[0 1.1],'Color',[0.8 0.8 0.8],'LineStyle','-'); end
    set(gca,'YTick',[0:0.2:1]); set(gca,'XTick',[-10 -5 0])
    box on
    ylabel('Frequency')
    xlabel('log2(H2BGFP int.)')

end
