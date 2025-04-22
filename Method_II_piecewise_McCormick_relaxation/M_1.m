%%% Manually convert the bilinear optimization problem into 2^4 = 16 by using Piecewise McCormick relaxation
%%% MATLAB toolbox YALMIP and the solver Mosek are used for calculating LPs 

%%% Setting the potential parameter
delta = ; 

%%% Set parameters and Calculate value of Bell ineuqality from the experimental data
t = ;
experiment_value = 0.0911 - t * (0.0015 + 0.0015 + 0.0008);
gamma = 0.5;

%%% Calculate value of Bell ineuqality for the \delta-OD model

%%% Define variables in the problem
a = sdpvar(1,2,'full'); % a(i+1):= p(a=0|x=i)
b = sdpvar(2,2,2,'full');% b(i+1,j+1,k+1):= p(b=0|x=i,y=j,a=k)
w = sdpvar(2,2,2,'full'); % w(i,j,k) := a(i) * b(i,j,k)

%%% Define objective function in the problem
% The objective funtion is the Bell inequality expression:
% p(a=b=0|x=y=1)-t*(p(a=1,b=0|x=0,y=1)+p(a=0,b=1|x=1,y=0)+p(a=b=0|x=y=0))
obj = w(2,2,1) - t * (a(2) - w(2,1,1) + b(1,2,2) - w(1,2,2) + w(1,1,1));

%%% Set conditions and Solve 16 LPs
opt_value = []; % Save the optimal value of 16 LPs in opt_value
for k = 0:15
    %%% Define conditions in the problem
    % None negative conditions for a
    cons=[0<=a(1)<=1;0<=a(2)<=1];
    % The PI (parameter independence) condition for Bob's output
    % \sum_a p(a,b=0|x=1,y) = \sum_a p(a,b=0|x=0,y), \forall y
    cons = [cons; w(2,2,1) + b(2,2,2) - w(2,2,2) == w(1,2,1) + b(1,2,2) - w(1,2,2)];
    cons = [cons; w(2,1,1) + b(2,1,2) - w(2,1,2) == w(1,1,1) + b(1,1,2) - w(1,1,2)];
    
    % Piecewise McCormick relaxation are applied on four variables p(b=0|x,y,a=1),i.e. b(x+1,y+1,2))
    % Each variable has two domain partitions [0,1/2] or [1/2,1]
    % Each combination of domain partitions is a bit string 0000-1111 (the MC_choice)
    MC_choice = dec2bin(k, 4) - '0';
    
    % For the current choice of MC_choice, set the correspounding conditions
    % for each relaxed variable b(x+1,y+1,2) and its related variables b(x+1,y+1,1)and w(x+1,y+1,a)
    for i = 1:4 
        index = dec2bin(i-1, 2)-'0';
        x = index(1); y = index(2);
        if MC_choice(i) == 0
            cons = [cons; 0 <= b(x+1,y+1,2) <= gamma];  
            cons = [cons; 0 <= b(x+1,y+1,1) <= gamma+delta]; 
            cons = [cons; condition_MC(x,y,0,w,a,b,gamma,delta)]; 
        elseif MC_choice(i) == 1
            cons = [cons; gamma <= b(x+1,y+1,2) <=1];        
            cons = [cons; gamma - delta <= b(x+1,y+1,1) <= 1];
            cons = [cons; condition_MC(x,y,1,w,a,b,gamma,delta)];
        end
    end
    
    %%% Solve the optimization problem
    ops = sdpsettings('solver', 'mosek'); 
    result = optimize(cons, - obj, ops); % YALMIP default optimization problem is minimization 
    
    %%% Save result for current choice of MC_choice
    if result.problem == 0
        opt_value(end + 1) = value(obj);  
    else
        yalmiperror(result.problem) % output error if any
        return
    end
end

% Output the biggest optimal value of the 16 LPs as the final optimal value 
final_opt_value = max(opt_value)

%%% Function to calculate the McCormick envelopes
function output = condition_MC(x,y,choice,w,a,b,gamma,delta)
La = 0; Ua = 1;
if choice == 0
    Lb = 0; Ub = gamma + delta;
    output1 = [w(x+1,y+1,1) >= Lb * a(x+1) + La * b(x+1,y+1,1) - La * Lb;
               w(x+1,y+1,1) <= Ub * a(x+1) + La * b(x+1,y+1,1) - La * Ub;
               w(x+1,y+1,1) <= Lb * a(x+1) + Ua * b(x+1,y+1,1) - Ua * Lb;
               w(x+1,y+1,1) >= Ub * a(x+1) + Ua * b(x+1,y+1,1) - Ua * Ub;];
    Lb = 0; Ub = gamma;
    output2 = [w(x+1,y+1,2) >= Lb * a(x+1) + La * b(x+1,y+1,2) - La * Lb;
               w(x+1,y+1,2) <= Ub * a(x+1) + La * b(x+1,y+1,2) - La * Ub;
               w(x+1,y+1,2) <= Lb * a(x+1) + Ua * b(x+1,y+1,2) - Ua * Lb;
               w(x+1,y+1,2) >= Ub * a(x+1) + Ua * b(x+1,y+1,2) - Ua * Ub;];
elseif choice == 1
    Lb = gamma - delta; Ub = 1;
    output1 = [w(x+1,y+1,1) >= Lb * a(x+1) + La * b(x+1,y+1,1) - La * Lb;
               w(x+1,y+1,1) <= Ub * a(x+1) + La * b(x+1,y+1,1) - La * Ub;
               w(x+1,y+1,1) <= Lb * a(x+1) + Ua * b(x+1,y+1,1) - Ua * Lb;
               w(x+1,y+1,1) >= Ub * a(x+1) + Ua * b(x+1,y+1,1) - Ua * Ub;];
    Lb = gamma; Ub = 1;
    output2 = [w(x+1,y+1,2) >= Lb * a(x+1) + Lb * b(x+1,y+1,2) - La * Lb;
               w(x+1,y+1,2) <= Ub * a(x+1) + La * b(x+1,y+1,2) - La * Ub;
               w(x+1,y+1,2) <= Lb * a(x+1) + Ua * b(x+1,y+1,2) - Ua * Lb;
               w(x+1,y+1,2) >= Ub * a(x+1) + Ua * b(x+1,y+1,2) - Ua * Ub;];
end
output = [output1; output2];
end