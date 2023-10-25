function fval = meg_modelComparison_rmse(sse1,sse2,k1,k2,n)
% F statistic to compare models based on sum of square errors 

%% Define variables
% 1 = simpler model 
% 2 = bigger model 
% sse2 will be smaller than sse1 with additional predictors 
% k2 = number of predictors (bigger model)
% k1 = number of predictors (simpler model) 
% sse2 = mse of bigger model 
% df = sample size - number of predictors 

%% 
n = 1; % sample size 
k2 = 4; 
k1 = 2; 
df1 = n - k1; 
df2 = n - k2; 

subjectIdx = 1; 
sse1 = mdlFit.linear.cueT1.session.fval(1,subjectIdx); 
sse2 = mdlFit.linear2Hz.cueT1.session.fval(1,subjectIdx); 

%% Calculate f statistic 
fval = ((sse1-sse2)/(k2-k1))/(sse2/(n-k2));

%% Calculate corresponding pval 
pval = 1-fcdf(fval,df1,df2); 

