
% load subject D2
% test criticality edge of chaos measure 
% 0 --> strong periodicity 
% 1 --> chaos 

% Estimate the chaoticity of a time-series y on a scale from
% zero to one, where zero indicates periodicity/stability and one indicates
% chaos/instability.

%% Create signal 
x = 1:0.0001:100; 
y = sin(x)*100; 


y = rand([1,1000]); % close to k 
%% Chaos test 
ds_method = 'minmax'; % 'minmax'
k = chaos_test(y,ds_method); 