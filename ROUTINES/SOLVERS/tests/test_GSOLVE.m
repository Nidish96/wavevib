clc
clear all
addpath('../')

x0 = [0;1];
%% With Fsolve
opts = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true);
Xs1 = fsolve(@(X) OBJFUN(X), x0, opts);
Xs1(2)

%% With GSOLVE (no scaling) 
opts = struct('Dscale', [1 1], 'Display', true, 'ITMAX', 100);
Xs2 = GSOLVE(@(X) OBJFUN(X), x0, opts);
Xs2(2)

%% With GSOLVE (with scaling & Line search)
opts = struct('Dscale', [1 1e-7], 'Display', true, 'ITMAX', 100, 'lsrch', 1, 'lsit', 30);
Xs3 = GSOLVE(@(X) OBJFUN(X), x0, opts);
Xs3(end)

%% Objective Function
function [R, dR] = OBJFUN(X)
    R = [(X(1)-1)^2;
        (1e7*X(2)-0.75)^2];
    dR = [2*(X(1)-1), 0;
        0, 2e7*(1e7*X(2)-0.75)];
end