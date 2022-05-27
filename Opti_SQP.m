%% EXAMPLE OF SQP ALGORITHM
clear variables; close all; clc;
fprintf("---------------------------------------------------------------\n")
fprintf("An implementation of Sequential Quadratic Programming method\nin a nonlinear constrained optimization problem\n")
fprintf("---------------------------------------------------------------\n")
%% INITIAL VALUES - INPUT
vars = 2; % number of variables
cons = 1; % number of constraints
maxIter=100; % max iterations
x = [-1;4]; % initial guess point
l =0; % LagrangeMultipliers vector
H = eye(vars,vars); % Hessian matrix assumed to be identity
%% EVALUATE AT STARTING POINT
fEval= f(x);
gEval = g(x);
[gViol,lViol] = Viols(gEval,l); 
gradfEval = gradf(x);
gradgEval = gradg(x);
P = Penalty(fEval,gViol,lViol);
%% SQP ALGORITHM
for i=1:maxIter
    %% SOLVE KKT CONDITIONS FOR THE OPTIMUM OF THE QUADRATIC APPROXIMATION
    sol = SolveKKT(gradfEval,gradgEval,gEval,H);
    xSol = sol(1:vars);
    lSol = sol(vars+1:vars+cons);
    %% IF THE LAGRANGE MULTIPLIER IS NEGATIVE SET IT TO ZERO
    for j = 1:length(lSol)
        if lSol(j)<0 
            sol= H\(-gradfEval)';
            xSol = sol(1:vars);
            lSol(j)=0;
        end
    end
    %% EVALUATE AT NEW CANDIDATE POINT
    xNew = x + xSol; 
    lNew = lSol;
    fEvalNew = f(xNew);
    gEvalNew = g(xNew);
    gradfEvalNew = gradf(xNew);
    gradgEvalNew = gradg(xNew);
    [gViolNew,lViolNew] = Viols(gEvalNew,lNew);
    PNew = Penalty(fEvalNew,gViolNew,lViolNew);
    %% IF PENALTY FUNCTION INCREASED, LOWER THE STEP BY 0.5
    while PNew-P>1e-4
        xSol = 0.5*xSol;
        xNew = x + xSol;
        fEvalNew = f(xNew);
        gEvalNew = g(xNew);
        gradfEvalNew = gradf(xNew);
        gradgEvalNew = gradg(xNew);
        [gViolNew,lViolNew] = Viols(gEvalNew,lNew);
        PNew = Penalty(fEvalNew,gViolNew,lViolNew);
    end
    %% STOPPING CRITERION
    if norm(xNew(1:vars)-x(1:vars))<=1e-2
        break
    end
    %% UPDATE THE HESSIAN
    gradLEval = gradLagr(gradfEval,gradgEval,lNew,vars); % lnew not l!!!
    gradLEvalNew = gradLagr(gradfEvalNew,gradgEvalNew,lNew,vars);
    Q = gradLEvalNew-gradLEval;
    dx = xNew-x;
    HNew = UpdateH(H,dx,Q);
    %% UPDATE ALL VALUES FOR NEXT ITERATION
    H = HNew;
    fEval = fEvalNew;
    gEval = gEvalNew;
    gradfEval = gradfEvalNew;
    gradgEval = gradgEvalNew;
    P = PNew;
    x = xNew;
end
fprintf('SQP: Optimum point:\n x1=%10.4f\n x2=%10.4f\n iterations =%10.0f \n', x(1), x(2), i)
%% FUNCTIONS NEEDED
function y = SolveKKT(gradfEval,gradgEval,gEval,Hessian)
A = [Hessian -gradgEval';gradgEval 0];
b = [-gradfEval -gEval]';
y = A\b;
end

function y = f(x)
y = x(1)^4 - 2*x(2)*x(1)^2 + x(2)^2 + x(1)^2 - 2*x(1)+5;
end

function y = gradf(x)
y(1) = 2*x(1)-4*x(1)*x(2)+4*x(1)^3-2;
y(2) = -2*x(1)^2 + 2*x(2);
end

function y = gradLagr(gradfEval,gradgEval,l,n)
y = gradfEval';
sum = zeros(n,1);
for i = 1:length(l)
    sum = sum -l(i)*gradgEval(i:n)';
end
y = y + sum;
end


function y = gradg(x)
y(1,1) = -2*x(1)-1/2;
y(1,2) = 3/4;
end

function y = g(x)
y(1) = -(x(1)+0.25)^2+0.75*x(2);
end

function [gViol,lViol] = Viols(gEval,l)
gViol=[];
lViol=[];
for i = 1:length(gEval)
    if gEval(i)<0
        gViol(i)=gEval(i);
        lViol(i)=l(i);
    end    
end
end


function y = Penalty(fEval,gViol,lViol)
sum = 0;
y = fEval;
for i = 1:length(gViol)
    sum = sum +  lViol(i)*abs(gViol(i));
end
y = y + sum;
end

function y = UpdateH(H,dx,gamma)
term1=(gamma*gamma') / (gamma'*dx);
term2 = ((H*dx)*(dx'*H)) / (dx'*(H*dx));
y = H + term1-term2;
end






