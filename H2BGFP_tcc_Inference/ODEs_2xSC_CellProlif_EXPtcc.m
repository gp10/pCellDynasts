function [f]=ODEs_2xSC_CellProlif_EXPtcc(t,x,lambdaS1,lambdaS2,gammaS2,mu)
%% Ordinary Differential Equations describing basal cell proliferation and suprabasal cell turnover according to 2xSC model (semi-coupled version)
% The frequency of basal and suprabasal cell pools experiencing different
% number of division rounds is resolved over time, as progenitor cells
% succesively divide in the basal layer and cells stratify and progress
% through suprabasal (spinous) compartment.
% This code implements default exponentially-distributed waiting time
% periods for all processes.

% from Piedrafita et al, 2020

%% Input:
% x: column vector containing the fraction of cell types found with each particular No. of division rounds.
% lambdaS1: average division rate of S1 population (/week)
% lambdaS2: average division rate of S2 population (/week)
% gammaS2: division-uncoupled stratification rate in S2 population (/week)
% mu: shedding rate from spinous layer (SL) or SL->GL transit time (/week)

%% Output:
% f: time-course changes in the fraction of cell types transitioning through different No. of division rounds.

%% Example:
% x = zeros(45,1); x([1,10,19,28,37]) = 1;
% lambdaS1 = 0.51*7; %/week
% lambdaS2 = 0.19*7; %/week
% gammaS2 = 0.19*7; %/week
% mu = 3.4596; %/week
% [t,x_div]=ode45(@(t,x) ODEs_2xSC_CellProlif_EXPtcc(t,x,lambdaS1,lambdaS2,gammaS2,mu),[0:20],x);

%% Initial definition of parameters:
f=zeros(27,1); 

%% Set of ODEs:
% f(1-9): S1div: 0-8
f(1) = -lambdaS1*x(1);
f(2) = lambdaS1*(x(1)-x(2));
f(3) = lambdaS1*(x(2)-x(3));
f(4) = lambdaS1*(x(3)-x(4));
f(5) = lambdaS1*(x(4)-x(5));
f(6) = lambdaS1*(x(5)-x(6));
f(7) = lambdaS1*(x(6)-x(7));
f(8) = lambdaS1*(x(7)-x(8));
f(9) = lambdaS1*(x(8)-x(9));
% f(10-18): S2div: 0-8
f(10) = - 0.5*(lambdaS2+gammaS2)*x(10);
f(11) = lambdaS2*x(10) - 0.5*(lambdaS2+gammaS2)*x(11);
f(12) = lambdaS2*x(11) - 0.5*(lambdaS2+gammaS2)*x(12);
f(13) = lambdaS2*x(12) - 0.5*(lambdaS2+gammaS2)*x(13);
f(14) = lambdaS2*x(13) - 0.5*(lambdaS2+gammaS2)*x(14);
f(15) = lambdaS2*x(14) - 0.5*(lambdaS2+gammaS2)*x(15);
f(16) = lambdaS2*x(15) - 0.5*(lambdaS2+gammaS2)*x(16);
f(17) = lambdaS2*x(16) - 0.5*(lambdaS2+gammaS2)*x(17);
f(18) = lambdaS2*x(17) - 0.5*(lambdaS2+gammaS2)*x(18);
% f(19-27): SLdiv: 0-8
f(19) = 0.5*gammaS2*x(10) - mu*x(19);
f(20) = lambdaS1*x(1) + 0.5*gammaS2*x(11) - mu*x(20);
f(21) = lambdaS1*x(2) + 0.5*gammaS2*x(12) - mu*x(21);
f(22) = lambdaS1*x(3) + 0.5*gammaS2*x(13) - mu*x(22);
f(23) = lambdaS1*x(4) + 0.5*gammaS2*x(14) - mu*x(23);
f(24) = lambdaS1*x(5) + 0.5*gammaS2*x(15) - mu*x(24);
f(25) = lambdaS1*x(6) + 0.5*gammaS2*x(16) - mu*x(25);
f(26) = lambdaS1*x(7) + 0.5*gammaS2*x(17) - mu*x(26);
f(27) = lambdaS1*x(8) + 0.5*gammaS2*x(18) - mu*x(27);
