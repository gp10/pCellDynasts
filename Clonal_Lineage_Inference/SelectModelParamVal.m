function [ParamVal] = SelectModelParamVal(selectModelInferMode,selectDataSet)
%% SELECT MODEL PARAMETER VALUES FROM A SPECIFIC SOURCE:
% A collection of parameter values obtained from different sources is
% presented. The function helps selecting the specific set of interest.

% from Piedrafita et al, 2020

%% Input:
% selectModelInferMode: string; source for parameter values, whether MLE (this work) or original publications
% selectDataSet: numeric; experimental data set the model parameter values should refer to

%% Output:
% ParamVal: structure containing the parameter values of interest

%%
if strcmpi('MLE',selectModelInferMode)

    switch selectDataSet

        case 1
            %% MLE PARAMETER VALUES FOR: OE_Lrig1_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 2.9; %(/week)
            ParamVal.r = 0.095;
            ParamVal.gamma = 5.3857; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 8;

        case 2
            %% MLE PARAMETER VALUES FOR: OE_Doupe2012_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 2.9; %(/week)
            ParamVal.r = 0.055;
            ParamVal.gamma = 3.6909; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 8;

        case 3
            %% MLE PARAMETER VALUES FOR: Paw_Lim2013_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 2.0; %(/week)
            ParamVal.r = 0.14;
            ParamVal.gamma = 2.2553; %(/week)
            ParamVal.tlag = 1.0/7; % one day
            ParamVal.GamShape = 4;

        case 4
            %% MLE PARAMETER VALUES FOR: Ear_Doupe2010_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 1.5; %(/week)
            ParamVal.r = 0.035;
            ParamVal.gamma = 1.7609; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 8;

        case 5
            %% MLE PARAMETER VALUES FOR: Back_Murai2018_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 1.2; %(/week)
            ParamVal.r = 0.04;
            ParamVal.gamma = 1.8769; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 16;

        case 6
            %% MLE PARAMETER VALUES FOR: Back_Fullgrabe2015_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 1.2; %(/week)
            ParamVal.r = 0.04;
            ParamVal.gamma = 1.8769; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 16;

        case 7
            %% MLE PARAMETER VALUES FOR: TailScale_SanchezDanes2016_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 1.2; %(/week)
            ParamVal.r = 0.14;
            ParamVal.gamma = 3.0857; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 8;

        case 8
            %% MLE PARAMETER VALUES FOR: TailInterscale_SanchezDanes2016_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 1.2; %(/week)
            ParamVal.r = 0.09;
            ParamVal.gamma = 2.2286; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 8;

    end

elseif strcmpi('Original',selectModelInferMode)

    switch selectDataSet

        case 1
            %% ORIGINAL PARAMETER VALUES FOR: OE_Lrig1_dataset
            disp('The parameter values reported are the MLE ones (data from this work)')
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 2.9; %(/week)
            ParamVal.r = 0.095;
            ParamVal.gamma = 5.3857; %(/week)
            ParamVal.tlag = 0.5/7; % half a day
            ParamVal.GamShape = 8;

        case 2
            %% ORIGINAL PARAMETER VALUES FOR: OE_Doupe2012_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 1.9; %(/week)
            ParamVal.r = 0.1;
            ParamVal.gamma = 3.5286; %(/week)
            ParamVal.tlag = 0/7;
            ParamVal.GamShape = 1; % EXP tcc

        case 3
            %% ORIGINAL PARAMETER VALUES FOR: Paw_Lim2013_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 2.0; %(/week)
            ParamVal.r = 0.2;
            ParamVal.gamma = 1.8; %(/week)
            ParamVal.tlag = 0/7;
            ParamVal.GamShape = 1; % EXP tcc

        case 4
            %% ORIGINAL PARAMETER VALUES FOR: Ear_Doupe2010_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 0.25; %(/week)
            ParamVal.r = 0.11;
            ParamVal.gamma = 0.0972; %(/week)
            ParamVal.tlag = 1/7; % one day
            ParamVal.GamShape = 1; % EXP tcc

        case 5
            %% ORIGINAL PARAMETER VALUES FOR: Back_Murai2018_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 1.16; %(/week)
            ParamVal.r = 0.06;
            ParamVal.gamma = 3.8835; %(/week)
            ParamVal.tlag = 0/7;
            ParamVal.GamShape = 1; % EXP tcc

        case 6
            %% ORIGINAL PARAMETER VALUES FOR: Back_Fullgrabe2015_dataset
            disp('No parameter fits were provided by Fullgrabe et al 2015')
            return

        case 7
            %% ORIGINAL PARAMETER VALUES FOR: TailScale_SanchezDanes2016_dataset
            ParamVal.modelType = 'SP'
            ParamVal.lambda = 1.2; %(/week)
            ParamVal.r = 0.22;
            ParamVal.gamma = 118.8; %(/week) (approximation for a value -> infinity)
            ParamVal.tlag = 0/7;
            ParamVal.GamShape = 1; % EXP tcc

        case 8
            %% ORIGINAL PARAMETER VALUES FOR: TailInterscale_SanchezDanes2016_dataset
            disp('Original fits in Sánchez-Danés et al 2016 using a SC-CP model')
            ParamVal.modelType = 'SC-CP'
            ParamVal.labelSfrac = 0.65;
            ParamVal.lambdaS = 0.45; %(/week)
            ParamVal.rS = 0.03;
            ParamVal.lambdaP = 1.7; %(/week)
            ParamVal.rP = 0.19;
            ParamVal.DeltaP = (ParamVal.rP+0.02) / ParamVal.rP - 1; %converts original formulation (from SD2016 Suppl.Methods, pg.3) to equivalent units for DeltaP as in Eq.4 of our Suppl.Methods
            ParamVal.gamma = 118.8; %(/week) (approximation for a value -> infinity)
            ParamVal.tlagS = 0/7;
            ParamVal.GamShapeS = 1; % EXP tcc
            ParamVal.tlagP = 0/7;
            ParamVal.GamShapeP = 1; % EXP tcc

    end
    
end
