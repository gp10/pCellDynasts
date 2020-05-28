function [LogLike_t, LogLike] = logLike_calc(rfreqs,myPDF,timepoints)
%% LOG-LIKELIHOOD CALCULATION OF THE MATCH PDF vs EMPIRICAL DISTRIBUTION:
% PDF corresponds with the numerical outcome of the probability density
% function for particular parameter values.

% from Piedrafita et al, 2020

%% Input:
% rfreqs: cell array of clone size frequencies, with format {:,timepoints}
% myPDF: matrix of clone size frequencies, with format (:,timepoints)
% timepoints: row-vector containing time points

%% Output:
% LogLike_t: 1xn vector of log-Likelihood values at the n individual time points
% LogLike: log-Likelihood value across time points

%%
LogLike_t = [];
LogLike = [];
for baz = 1:length(timepoints)
    if (size(myPDF,1) < size(rfreqs{:,baz},1))
        scale_up = size(rfreqs{:,baz},1);
        myPDF = [myPDF; zeros(scale_up-size(myPDF,1),size(myPDF,2))];
    end
    % Calculation of log-Likelihood value for each time point:
    LogLike_t(1,baz) = sum ( rfreqs{:,baz}(find(rfreqs{:,baz}~=0),1) .* log(myPDF(find(rfreqs{:,baz}~=0),baz)) );
end
% Calculation of log-Likelihood value across time points:
LogLike = sum(LogLike_t,2);
