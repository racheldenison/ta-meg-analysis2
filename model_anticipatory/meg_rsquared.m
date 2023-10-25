function rsq = meg_rsquared(ssr,sst)
% Calculate r-squared 
% Inputs:
%   ssr: sum squared regression, scaler 
%   sst: total sum of squares, (number of samples), scaler
% Ouputs 
%   rsq: r-squared, scaler 

rsq = 1 - ssr/sst; 