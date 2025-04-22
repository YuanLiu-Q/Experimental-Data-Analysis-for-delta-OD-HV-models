%%% Solve the bilinear optimization problem by using Gurobi Optimizer
%%% Gurobi Optimizer is interfaced by the MATLAB toolbox YALMIP

%%% Setting the potential parameter
delta = ;

%%% Load experimental data and calculate value of Bell ineuqality from the data
M = ;
t = ;
experiment_value = Cal_Exp_Data(M,t); 

%%% Calculate value of Bell ineuqality for the \delta-OD model

%%% Define variables in the problem
N = M + 1; %The number of measurement settings for each party is M+1
a = sdpvar(1,N,'full'); % a(i+1):= p(a=0|x=i)
b = sdpvar(N,N,2,'full'); % b(i+1,j+1,k+1):= p(b=0|x=i,y=j,a=k)

%%% Define objective function in the problem
% The objective funtion is the Bell inequality expression:
% p(a=b=0|x=y=M)-t*\sum_{i=0}^{M-1}(p(a=1,b=0|x=i,y=i+1)+p(a=0,b=1|x=i+1,y=i))-t*p(a=b=0|x=y=0)
obj = a(N) * b(N,N,1) - t * a(1) * b(1,1,1);
for i = 0 : M-1
    obj = obj - t * (1 - a(i+1)) * b(i+1,i+2,2) - t * a(i+2) * (1 - b(i+2,i+1,1));
end

%%% Define conditions in the problem
cons = [];
% None negative conditions
cons = [cons; 0 <= a <= 1; 0 <= b <= 1];
% The delta-OD (output-dependence) conditions
for i = 1 : N
    for j = 1 : N
        cons = [cons; b(i,j,2) - delta <= b(i,j,1) <= b(i,j,2) + delta];
    end
end
% The PI (parameter independence) condition for Bob's output
for y = 1 : N % change index of y
    for x = 2 : N % the PI condition for a given y
    % \sum_a p(a,b=0|x,y) = \sum_a p(a,b=0|x=0,y), \forall x,y
    cons = [cons; a(x) * b(x,y,1) + (1 - a(x)) * b(x,y,2) == a(1) * b(1,y,1) + (1 - a(1)) * b(1,y,2)];
    end
end

%%% Solve the optimization problem
ops = sdpsettings('solver','gurobi', ...
    'gurobi.Cuts', 2, ...
    'gurobi.MIPFocus', 3, ...
    'gurobi.CutPasses', 8, ...
    'gurobi.GomoryPasses', 3, ...
    'gurobi.CliqueCuts', 2, ...
    'gurobi.Presolve', 0, ...
    'gurobi.MIPGap', 0.05, ...
    'gurobi.Heuristics', 0.1,...
    'gurobi.TuneTimeLimit',0,...
    'gurobi.Method', 2, ...
    'gurobi.Crossover', 0);
result = optimize(cons, - obj, ops); % YALMIP default optimization problem is minimization 

%%% Output result
if result.problem == 0
    value(obj) % output the optimal value of objective function
else
    yalmiperror(result.problem) % output error if any
    
end



%%% Function to calculate value of Bell ineuqality from the data
function output=Cal_Exp_Data(M,t)
    if M==20
        output=0.275;
    elseif M==25
        output=0.257;
    elseif M==30
        output=0.199;
    else
        if M == 1
        data = readtable('Data_Stockholm.xlsx');
        elseif M == 5
        data = readtable('Data_Oxford.xlsx');
        elseif M == 40
            data = readtable('Data_Roma_A.xlsx');
        elseif M == 44
            data = readtable('Data_Roma_B.xlsx');
        else
            disp('Missing data for M');
            return
        end
    col1 = data{:,1};
    col1 = col1(~isnan(col1));
    col2 = data{:,2};
    col2 = col2(~isnan(col2));
    output = sum(col1) - t * sum(col2);
    end
end