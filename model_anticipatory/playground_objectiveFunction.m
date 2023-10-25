function [E, yhat1, yhat2] = playground_objectiveFunction(x,data1,data2,t,paramNames)

%% 
idx = contains(paramNames,'amplitude1');
amplitude1 = x(idx);
idx = contains(paramNames,'amplitude2');
amplitude2 = x(idx);

%% 
yhat1 = amplitude1 * sin(t);
yhat2 = amplitude2 * cos(t);

y1 = data1;
y2 = data2;

%% 
E = sum ( (yhat1-y1).^2 + (yhat2-y2).^2 );
nVars = numel(paramNames); 

