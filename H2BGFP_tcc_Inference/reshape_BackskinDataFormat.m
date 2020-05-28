function [new_rtime,new_nmice,new_H2BGFP,new_FieldView,new_log2_H2BGFP,new_nlog2_H2BGFP,new_mean_log2_H2BGFP,new_mean_nlog2_H2BGFP,new_H2BGFP_all,new_nlog2_H2BGFP_all] = reshape_BackskinDataFormat(rtime,nmice,H2BGFP,FieldView,log2_H2BGFP,nlog2_H2BGFP,mean_log2_H2BGFP,mean_nlog2_H2BGFP,H2BGFP_all,nlog2_H2BGFP_all)
%% Algorithm used to reshape back-skin H2BGFP original variables according to time sorting
% Back skin H2BGFP data come from 3 different experiments and were
% originally organized accordingly, separately for the 3 sets. This script
% reshapes H2BGFP-related variables so as to organize data according to
% collection time, merging results from the different experiments whenever
% is appropriate.

% from Piedrafita et al, 2020

%% Input:
% rtime: horizontal vector containing experimental time points (expressed in weeks); contains repeats as they stand for 3 different experiments
% *: other variables whose format is to be changed according to map time points in the given order (variable type is preserved)

%% Output:
% new_rtime: new vector of time points after merging data from different experiments collected at the same time and sorting other time points
% new_*: other variable formats are reshaped according with the same time pattern (preserving their variable type).

%% Initial definition of parameters:
new_H2BGFP = {}; new_FieldView = {}; new_log2_H2BGFP = {}; new_nlog2_H2BGFP = {};
new_mean_log2_H2BGFP = {}; new_mean_nlog2_H2BGFP = {}; new_H2BGFP_all = {}; new_nlog2_H2BGFP_all = {};
new_nmice = 0;
new_rtime = [];

%% UPDATE TIME-POINT VECTOR:
new_rtime = rtime([1,2,5,3,7]);

%% MERGE DATA FROM t=0
counter = 0;
new_mean_log2_H2BGFP_t0 = []; new_mean_nlog2_H2BGFP_t0 = []; new_H2BGFP_all_t0 = []; new_nlog2_H2BGFP_all_t0 = [];
for aja = [1,4,6] % data elements 1,4,6 in input variables go into position 1 in output variables
    for eje = 1:nmice(aja)
        counter = counter + 1;
        new_H2BGFP{counter,1} = H2BGFP{eje,aja}; new_FieldView{counter,1} = FieldView{eje,aja};
        new_log2_H2BGFP{counter,1} = log2_H2BGFP{eje,aja};
        new_nlog2_H2BGFP{counter,1} = nlog2_H2BGFP{eje,aja};
    end
    new_mean_log2_H2BGFP_t0 = [new_mean_log2_H2BGFP_t0 mean_log2_H2BGFP{1,aja}];
    new_mean_nlog2_H2BGFP_t0 = [new_mean_nlog2_H2BGFP_t0 mean_nlog2_H2BGFP{1,aja}];
    new_H2BGFP_all_t0 = [new_H2BGFP_all_t0; H2BGFP_all{1,aja}];
    new_nlog2_H2BGFP_all_t0 = [new_nlog2_H2BGFP_all_t0; nlog2_H2BGFP_all{1,aja}];
    new_nmice = new_nmice + nmice(aja);
end
new_mean_log2_H2BGFP{1,1} = new_mean_log2_H2BGFP_t0; new_mean_nlog2_H2BGFP{1,1} = new_mean_nlog2_H2BGFP_t0; new_H2BGFP_all{1,1} = new_H2BGFP_all_t0; new_nlog2_H2BGFP_all{1,1} = new_nlog2_H2BGFP_all_t0;

%% SORT OTHER DATA ELEMENTS BY TIME:
counti = 0;
for aja = [2,5,3,7] % data elements 2,5,3,7 in input variables are sorted in output variables according to acquisition time
    counti = counti + 1;
    for eje = 1:nmice(aja)
        new_H2BGFP{eje,counti+1} = H2BGFP{eje,aja}; new_FieldView{eje,counti+1} = FieldView{eje,aja};
        new_log2_H2BGFP{eje,counti+1} = log2_H2BGFP{eje,aja};
        new_nlog2_H2BGFP{eje,counti+1} = nlog2_H2BGFP{eje,aja};
    end
    new_mean_log2_H2BGFP{1,counti+1} = mean_log2_H2BGFP{1,aja};
    new_mean_nlog2_H2BGFP{1,counti+1} = mean_nlog2_H2BGFP{1,aja};
    new_H2BGFP_all{1,counti+1} = H2BGFP_all{1,aja};
    new_nlog2_H2BGFP_all{1,counti+1} = nlog2_H2BGFP_all{1,aja};
    new_nmice(1,counti+1) = nmice(1,aja);
end
