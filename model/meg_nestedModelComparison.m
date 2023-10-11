function fval = meg_nestedModelComparison(sse1,sse2,k1,k2,n)
% Calculate F-statistic to compare nested models 
% Inputs:
%   see: sum of squared errors (scaler) 
%       -- sse1: for simpler model 
%       -- sse2: for bigger model (decrease in error with additional
%       predictors) 
%   k: number of parameters (scaler) 
%       -- k1: for simpler model
%       -- k2: for bigger model 
%   n: sample size, # of time points 
% Ouputs 
%   fval: F-statistic
%   -- Does adding additional predictors significantly improve the fit of the model? 

%% 
fval = ( (sse1-sse2)/(k2-k1) ) / ( sse2/(n - k2) );