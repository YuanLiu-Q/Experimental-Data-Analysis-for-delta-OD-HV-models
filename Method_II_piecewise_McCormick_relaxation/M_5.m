%%% Manually convert the bilinear optimization problem into 2^(12) = 4096 by using Piecewise McCormick relaxation
%%% MATLAB toolbox YALMIP and the solver Mosek are used for calculating LPs 

%%% Setting the potential parameter
delta = ;  

%%% Set parameters and Calculate value of Bell ineuqality from the experimental data
t = ;
experiment_value = 0.55 - t * (0.051 + 0.052 + 0.058 + 0.029 + ...
                   0.034 + 0.029 + 0.045 + 0.017 + 0.027 + 0.014 + 0.024);
gamma = 0.5;

%%% Calculate value of Bell ineuqality for the \delta-OD model

%%% Define variables in the problem
a = sdpvar(1,6,'full'); % a(i+1):= p(a=0|x=i)
b = sdpvar(6,6,2,'full');% b(i+1,j+1,k+1):= p(b=0|x=i,y=j,a=k)
w = sdpvar(6,6,2,'full'); % w(i,j,k) := a(i) * b(i,j,k)

%%% Define objective function in the problem
% The objective funtion is the Bell inequality expression:
% p(a=b=0|x=y=5)-t*\sum_{i=0}^{4}(p(a=1,b=0|x=i,y=i+1)+p(a=0,b=1|x=i+1,y=i))-t*p(a=b=0|x=y=0)
obj = w(6,6,1) - t * w(1,1,1);
for i = 0:4
    obj = obj - t * (b(i+1,i+2,2) - w(i+1,i+2,2)) - t * (a(i+2) - w(i+2,i+1,1));
end

%%% Set conditions and Solve 4096 LPs

% Piecewise McCormick relaxation are only applied on 12 variables
% p(b=0|x,y,a=1),i.e. b(x+1,y+1,2)), for certain (x,y) pairs listed here
Relax_xy = [0,0; 0,1; 1,0; 1,2; 2,1; 2,3; 3,2; 3,4; 4,3; 4,5; 5,4; 5,5;];
opt_value = []; % Save the optimal value of 4096 LPs in opt_value
for k = 0 : 2^(12)-1
    %%% Define conditions in the problem
    % None negative conditions for a
    cons = [0 <= a <= 1];
    % Conditions for non-relaxed variables b 
    for i = 1:6
        for j = 1:6
            if isempty(find(all(Relax_xy == [i,j],2)))
                % None negative conditions for b
                cons = [cons; 0 <= b(i,j,1) <= 1; 0 <= b(i,j,2) <=1];
                % The delta-OD (output-dependence) conditions for b
                cons = [cons; - delta + b(i,j,2) <= b(i,j,1) <= delta + b(i,j,2)];
            end
        end
    end
    % The PI (parameter independence) condition for Bob's output
    for y = 1:6 % change index of y
        for x = 2:6 % the PI condition for a given y
        % \sum_a p(a,b=0|x,y) = \sum_a p(a,b=0|x=0,y), \forall x,y
        cons = [cons; w(x,y,1) + b(x,y,2) - w(x,y,2) == w(1,y,1) + b(1,y,2) - w(1,y,2)];
        end
    end
    % Piecewise McCormick relaxation are applied on 12 variables p(b=0|x,y,a=1),i.e. b(x+1,y+1,2)
    % With certain (x,y) pairs listed in Relax_xy
    % Each variable has two domain partitions [0,1/2] or [1/2,1]
    % Each combination of domain partitions is a bit string of length 12 (the MC_choice)
    MC_choice = dec2bin(k, 12) - '0'; 
    % For the current choice of MC_choice, set the correspounding conditions
    % for each relaxed variable b(x+1,y+1,2) and its related variables b(x+1,y+1,1) and w(x+1,y+1,a)
    for i = 1:12
        index = Relax_xy(i,:);
        x = index(1); y = index(2);
        if MC_choice(i) == 0
            cons = [cons; 0 <= b(x+1,y+1,2) <= gamma];  
            cons = [cons; 0 <= b(x+1,y+1,1) <= gamma+delta]; 
            cons = [cons; condition_MC(x,y,0,w,a,b,gamma,delta)];
        elseif MC_choice(i) == 1
            cons = [cons; gamma <= b(x+1,y+1,2) <= 1];        
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

% Output the biggest optimal value of the 4096 LPs as the final optimal value
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
    Lb = gamma-delta; Ub = 1;
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
