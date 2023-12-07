function [rsq rsq2 y_fit_meanCorrected] = calculateRSQ(y,y_fit,meanSub)
% function rsq = calculateRSQ(y,y_fit,meanSub)
% y: real data 
% y_fit: fitted data 
% Calculate r-squared 

if meanSub % mean subtraction 
    meanDiff = mean(y,'omitnan')-mean(y_fit,'omitnan'); 
    % err = y-y_fit+meanDiff;  
    y_fit_meanCorrected = y_fit+meanDiff; 
    err = y_fit_meanCorrected - y;  
else
    err=y-y_fit;
end

SSres=sum(err.^2);
SStot=sum((y-mean(y)).^2);
rsq=1-(SSres./SStot);

%% alt
r = corrcoef(y,y_fit);
rsq2 = r(2)^2; 

