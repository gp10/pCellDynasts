function [f]=ODEs_2xSC_CellProlif_GAMtcc(t,x,lambdaS1,lambdaS2,gammaS2,mu)
%% Ordinary Differential Equations describing basal cell proliferation and suprabasal cell turnover according to 2xSC model (semi-coupled version)
% The frequency of basal and suprabasal cell pools experiencing different
% number of division rounds is resolved over time, as progenitor cells
% succesively divide in the basal layer and cells stratify and progress
% through suprabasal (spinous) compartment.
% This code implements gamma-distributed waiting time periods for the S1
% population division process (cell cycle time) and SL shedding process
% (suprabasal transit time), both with Gamma shape parameter kappa=2. This
% is done by replacing these Poisson processes by two-step processes
% involving two subpopulations, each showing half the mean waiting time.

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
% [t,x_div]=ode45(@(t,x) ODEs_2xSC_CellProlif_GAMtcc(t,x,lambdaS1,lambdaS2,gammaS2,mu),[0:20],x);

%% Initial definition of parameters:
f=zeros(45,1);

%% Set of ODEs:
% f(1-9): S1a_div: 0-8
f(1) = -2*lambdaS1*x(1);
f(2) = 2*lambdaS1*(x(10)-x(2));
f(3) = 2*lambdaS1*(x(11)-x(3));
f(4) = 2*lambdaS1*(x(12)-x(4));
f(5) = 2*lambdaS1*(x(13)-x(5));
f(6) = 2*lambdaS1*(x(14)-x(6));
f(7) = 2*lambdaS1*(x(15)-x(7));
f(8) = 2*lambdaS1*(x(16)-x(8));
f(9) = 2*lambdaS1*(x(17)-x(9));
% f(10-18): S1b_div: 0-8
f(10) = 2*lambdaS1*(x(1)-x(10));
f(11) = 2*lambdaS1*(x(2)-x(11));
f(12) = 2*lambdaS1*(x(3)-x(12));
f(13) = 2*lambdaS1*(x(4)-x(13));
f(14) = 2*lambdaS1*(x(5)-x(14));
f(15) = 2*lambdaS1*(x(6)-x(15));
f(16) = 2*lambdaS1*(x(7)-x(16));
f(17) = 2*lambdaS1*(x(8)-x(17));
f(18) = 2*lambdaS1*(x(9)-x(18));
% f(19-27): S2_div: 0-8
f(19) = - 0.5*(lambdaS2+gammaS2)*x(19);
f(20) = lambdaS2*x(19) - 0.5*(lambdaS2+gammaS2)*x(20);
f(21) = lambdaS2*x(20) - 0.5*(lambdaS2+gammaS2)*x(21);
f(22) = lambdaS2*x(21) - 0.5*(lambdaS2+gammaS2)*x(22);
f(23) = lambdaS2*x(22) - 0.5*(lambdaS2+gammaS2)*x(23);
f(24) = lambdaS2*x(23) - 0.5*(lambdaS2+gammaS2)*x(24);
f(25) = lambdaS2*x(24) - 0.5*(lambdaS2+gammaS2)*x(25);
f(26) = lambdaS2*x(25) - 0.5*(lambdaS2+gammaS2)*x(26);
f(27) = lambdaS2*x(26) - 0.5*(lambdaS2+gammaS2)*x(27);
% f(28-36): SLa_div: 0-8
f(28) = 0.5*gammaS2*x(19) - 2*mu*x(28);
f(29) = 2*lambdaS1*x(10) + 0.5*gammaS2*x(20) - 2*mu*x(29);
f(30) = 2*lambdaS1*x(11) + 0.5*gammaS2*x(21) - 2*mu*x(30);
f(31) = 2*lambdaS1*x(12) + 0.5*gammaS2*x(22) - 2*mu*x(31);
f(32) = 2*lambdaS1*x(13) + 0.5*gammaS2*x(23) - 2*mu*x(32);
f(33) = 2*lambdaS1*x(14) + 0.5*gammaS2*x(24) - 2*mu*x(33);
f(34) = 2*lambdaS1*x(15) + 0.5*gammaS2*x(25) - 2*mu*x(34);
f(35) = 2*lambdaS1*x(16) + 0.5*gammaS2*x(26) - 2*mu*x(35);
f(36) = 2*lambdaS1*x(17) + 0.5*gammaS2*x(27) - 2*mu*x(36);
% f(37-45): SLb_div: 0-8
f(37) = 2*mu*(x(28)-x(37));
f(38) = 2*mu*(x(29)-x(38));
f(39) = 2*mu*(x(30)-x(39));
f(40) = 2*mu*(x(31)-x(40));
f(41) = 2*mu*(x(32)-x(41));
f(42) = 2*mu*(x(33)-x(42));
f(43) = 2*mu*(x(34)-x(43));
f(44) = 2*mu*(x(35)-x(44));
f(45) = 2*mu*(x(36)-x(45));
