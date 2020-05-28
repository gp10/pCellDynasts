function [ParamVal] = SelectModelParamVal(selectModelInferMode,selectDataSet)
%% SELECT MODEL PARAMETER VALUES FROM A SPECIFIC SOURCE:
% A collection of parameter values obtained from different sources is
% presented. The function helps selecting the specific set of interest.

% from Piedrafita et al, 2020

%% Input:
% selectModelInferMode: string; source for parameter values, whether ABC (inference methods in this work) or original publications/ad-hoc
% selectDataSet: numeric; experimental data set or condition model parameter values map to

%% Output:
% ParamVal: structure containing the model/parameter values of interest

%%
if strcmpi('ABC_inference',selectModelInferMode)

    switch selectDataSet

        case 1
            %% OPTIMAL tcc PARAMETER VALUES FOR: OE_R26rtTA
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 2.9; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 8;

        case 2
            %% OPTIMAL tcc PARAMETER VALUES FOR: Paw_R26rtTA
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 2.0; %(/week)
            ParamVal.tlag = 1.0/7; % one day
            ParamVal.GamShape = 4;

        case 3
            %% OPTIMAL tcc PARAMETER VALUES FOR: Ear_R26rtTA
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 1.5; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 8;

        case 4
            %% OPTIMAL tcc PARAMETER VALUES FOR: Back_R26rtTA
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 1.2; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 16;
            
        case 5
            %% OPTIMAL tcc PARAMETER VALUES FOR: Tail_Krt5_tTA_Mascre2012
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 1.337; %(/week)
            ParamVal.tlag = 0.6/7; % a bit more than half a day
            ParamVal.GamShape = 2^[0.5];
            
        case 6
            %% OPTIMAL tcc PARAMETER VALUES FOR: Back_Krt5_tTA_Sada2016
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 3.1104; %(/week)
            ParamVal.tlag = 0.4/7; % a bit less than half a day
            ParamVal.GamShape = 1;
            ParamVal.mu = 49.7665; %(/week)
            ParamVal.tlagShed = 0;
            ParamVal.GamShapeShed = 2^[1.5];
            
    end

elseif strcmpi('Original/Ad-hoc',selectModelInferMode)

    switch selectDataSet

        case 1
            %% AD-HOC SP MODEL WITH SOME REALISTIC PARAMETER VALUES
            disp('Example of realistic SP model')
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 2.0; %(/week)
            ParamVal.r = 0.1;
            ParamVal.gamma = 3.0; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 8;

        case 2
            %% ORIGINAL PARAMETER VALUES IN Mascre et al: Tail_Mascre2012_dataset (with reasonable assumptions on tcc distribution shape)
            disp('SC-CP model from Mascré et al 2012')
            ParamVal.modelType = 'SC-CP'
            ParamVal.lambdaS = 0.096; %(/week) (= 5/year)
            ParamVal.rS = 0.1;
            ParamVal.lambdaP = 1.21; %(/week)
            ParamVal.rP = 0.1;
            ParamVal.DeltaP = 0.14/2; %converts original formulation (from Mascre2012 Suppl.Methods, pg.18) to equivalent units for DeltaP as in Eq.4 of our Suppl.Methods
            ParamVal.gamma = 4.8; %(/week)
            omega = 2*ParamVal.lambdaP*ParamVal.DeltaP*ParamVal.rP;
            ParamVal.labelSfrac = omega / (omega + ParamVal.lambdaS); % see Eq.7 of our Suppl.Methods (assuming D-cells are not induced)
            ParamVal.tlagS = 7/7; % reasonable assumption
            ParamVal.GamShapeS = 8; % GAM tcc
            ParamVal.tlagP = 0.5/7; % reasonable assumption
            ParamVal.GamShapeP = 8; % GAM tcc

        case 3
            %% ORIGINAL PARAMETER VALUES IN SanchezDanes et al: TailInterscale_SanchezDanes2016_dataset (with reasonable assumptions on tcc distribution shape)
            disp('SC-CP model from Sánchez-Danés et al 2016')
            ParamVal.modelType = 'SC-CP'
            ParamVal.lambdaS = 0.45; %(/week)
            ParamVal.rS = 0.03;
            ParamVal.lambdaP = 1.7; %(/week)
            ParamVal.rP = 0.19;
            ParamVal.DeltaP = (ParamVal.rP+0.02) / ParamVal.rP - 1; %converts original formulation (from SD2016 Suppl.Methods, pg.3) to equivalent units for DeltaP as in Eq.4 of our Suppl.Methods
            ParamVal.gamma = 300; %(/week) (approximation for a value -> infinity, as assumed in SD2016)
            omega = 2*ParamVal.lambdaP*ParamVal.DeltaP*ParamVal.rP;
            ParamVal.labelSfrac = omega / (omega + ParamVal.lambdaS); % see Eq.7 of our Suppl.Methods (assuming D-cells are not induced)
            ParamVal.tlagS = 2/7; % reasonable assumption
            ParamVal.GamShapeS = 8; % GAM tcc
            ParamVal.tlagP = 0.5/7; % reasonable assumption
            ParamVal.GamShapeP = 8; % GAM tcc

        case 4
            %% ORIGINAL PARAMETER VALUES IN Sada et al: Back_Sada2016_dataset (with reasonable assumptions on tcc distribution shape)
            disp('2xSC model (hybrid version) from Sada et al 2016')
            ParamVal.modelType = '2xSC'
            ParamVal.lambdaS1 = 0.47*7; %/week
            ParamVal.lambdaS2 = 0.19*7; %/week
            ParamVal.uS1 = 0.20;
            ParamVal.gammaS1 = ParamVal.lambdaS1; % to satisfy homeostasis
            ParamVal.gammaS2 = ParamVal.lambdaS2; % to satisfy homeostasis
            ParamVal.labelS1frac = 0.74;
            ParamVal.densSL_vs_densBL = 0.78; % experimentally determined
            ParamVal.mu = (ParamVal.lambdaS1*ParamVal.labelS1frac/(1+ParamVal.uS1) + ParamVal.lambdaS2*(1-ParamVal.labelS1frac)*0.5)/ParamVal.densSL_vs_densBL; % to satisfy homeostasis (pg.32 of our Suppl.Methods)
            ParamVal.tlagDivS1 = 0.5/7; % reasonable assumption
            ParamVal.GamShapeDivS1 = 8; % GAM division time
            ParamVal.tlagStrS1 = 0;
            ParamVal.GamShapeStrS1 = 1; % EXP stratif time
            ParamVal.tlagDivS2 = 0.5/7; % reasonable assumption
            ParamVal.GamShapeDivS2 = 8; % GAM division time
            ParamVal.tlagStrS2 = 0;
            ParamVal.GamShapeStrS2 = 1; % EXP stratif time
            ParamVal.tlagShed = 0;
            ParamVal.GamShapeShed = 1; % EXP shedding time
            
        case 5
            %% ORIGINAL PARAMETER VALUES IN Sada et al: Back_Sada2016_dataset
            disp('2xSC model (semi-coupled version) from Sada et al 2016 parameters')
            ParamVal.modelType = '2xSC'
            ParamVal.lambdaS1 = 0.51*7; %/week
            ParamVal.lambdaS2 = 0.19*7; %/week
            ParamVal.uS1 = 0;
            ParamVal.gammaS1 = ParamVal.lambdaS1; % to satisfy homeostasis
            ParamVal.gammaS2 = ParamVal.lambdaS2; % to satisfy homeostasis
            ParamVal.labelS1frac = 0.70;
            ParamVal.densSL_vs_densBL = 0.78; % experimentally determined
            ParamVal.mu = (ParamVal.lambdaS1*ParamVal.labelS1frac/(1+ParamVal.uS1) + ParamVal.lambdaS2*(1-ParamVal.labelS1frac)*0.5)/ParamVal.densSL_vs_densBL; % to satisfy homeostasis (pg.32 of our Suppl.Methods)
            ParamVal.tlagDivS1 = 0;
            ParamVal.GamShapeDivS1 = 2;
            ParamVal.tlagStrS1 = 0;
            ParamVal.GamShapeStrS1 = 2;
            ParamVal.tlagDivS2 = 0;
            ParamVal.GamShapeDivS2 = 1;
            ParamVal.tlagStrS2 = 0;
            ParamVal.GamShapeStrS2 = 1;
            ParamVal.tlagShed = 0;
            ParamVal.GamShapeShed = 2;
            
        case 6
            %% PARAMETER VALUES INFERRED FOR Sada et al's BEST 2xSC MODEL FIT: Back_Sada2016_dataset (pg.33 of our Suppl.Methods)
            % (these parameter values account for the probability terms ignored in Sada et al's equations)
            disp('2xSC model (semi-coupled version) from Sada et al 2016')
            ParamVal.modelType = '2xSC'
            ParamVal.lambdaS1 = 0.51*7; %/week
            ParamVal.lambdaS2 = 0.38*7; %/week
            ParamVal.uS1 = 0;
            ParamVal.gammaS1 = ParamVal.lambdaS1; % to satisfy homeostasis
            ParamVal.gammaS2 = ParamVal.lambdaS2; % to satisfy homeostasis
            ParamVal.labelS1frac = 0.70;
            ParamVal.densSL_vs_densBL = 0.78; % experimentally determined
            ParamVal.mu = (ParamVal.lambdaS1*ParamVal.labelS1frac + ParamVal.lambdaS2*(1-ParamVal.labelS1frac))/ParamVal.densSL_vs_densBL; % miscalculating homeostasis requirement (see pg.32-33 of our Suppl.Methods)
            ParamVal.tlagDivS1 = 0;
            ParamVal.GamShapeDivS1 = 2;
            ParamVal.tlagStrS1 = 0;
            ParamVal.GamShapeStrS1 = 2;
            ParamVal.tlagDivS2 = 0;
            ParamVal.GamShapeDivS2 = 1;
            ParamVal.tlagStrS2 = 0;
            ParamVal.GamShapeStrS2 = 1;
            ParamVal.tlagShed = 0;
            ParamVal.GamShapeShed = 2;
            
    end
    
end
