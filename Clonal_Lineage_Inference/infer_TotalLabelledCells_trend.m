function [avgNCells, stdNCells] = infer_TotalLabelledCells_trend(avgCloneSize,semCloneSize,avgNClones,stdNClones)
%% INFERENCE OF PERCENTAGE OF LABELLED CELLS BASED ON AVERAGE CLONE SIZE AND CLONE DENSITY VALUES:
% The percentage of labelled basal cells is inferred as a product of the
% average basal clone size multiplied by the clone density or total number
% of clones. However, as these two variables are measured from independent
% experiments, error propagation has to be performed to account for the
% level of uncertainty in the inferred value of %labelled cells.

% R = A * B      | A=mean(CloneSize) ; B=NClones
% dR = dA*B + dB*A + dA*dB      | dA=sem(CloneSize) ; dB=std(NClones)

% from Piedrafita et al, 2020

%% Input:
% avgCloneSize: row vector containing average clone size values at different time points
% semCloneSize: standard error of the mean clone size for same time points
% avgNClones: row vector containing average number of clones per unit area for each time point
% stdNClones: standard deviation of the number of clones per unit area for each time point

%% Output:
% avgNCells: inferred number of labelled basal cells per unit area for the specified time points
% stdNCells: inferred standard deviation of the number of labelled basal cells per unit area for the specified time points

%%
% Calculation of avgNCells:
avgNCells = avgCloneSize .* avgNClones;
% Calculation of stdNCells (by error propagation, as described above):
stdNCells = semCloneSize .* avgNClones + stdNClones .* avgCloneSize + semCloneSize .* stdNClones;
