function rsq = calculateRSQ(y,y_fit,meanSub)
% function rsq = calculateRSQ(y,y_fit,meanSub)
% y: real data 
% y_fit: fitted data 
% Calculate r-squared 

if meanSub % mean subtraction 
    meanDiff = mean(y_fit,'omitnan')-mean(y,'omitnan'); 
    err = y-y_fit+meanDiff;  
else
    err=y-y_fit;
end

SSres=sum(err.^2);
SStot=sum((y-mean(y)).^2);
rsq=1-(SSres./SStot);