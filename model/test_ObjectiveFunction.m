

t = 0:1:1000; 

Fs = 1000; 

data = 5* sin((2*pi/(500))*t+300); 
noise = -0.5 + rand(size(data)); % uniformly distributed random noise 
data = data + noise; 

% [E] = objectiveFunction(x,data,t,Fs); 

x0(1) = 0; % starting coefficient 
x0(2) = 0; 

lb(1) = 0; 
ub(2) = Inf; 

% type objectiveFunction
fun = @(x)objectiveFunction(x,data,t);
[solution,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);

figure
hold on 
plot(t,data,'LineWidth',3)
[E,yhat] = objectiveFunction(solution,data,t,Fs); 
plot(t,yhat,'LineWidth',1)