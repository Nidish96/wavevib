clc
clear all
addpath('../')

%% 
mul = 100;
po = 2;

c1 = 4*mul;
c2 = 1*mul;

%% 
xstart = -5*mul;
xend = 5*mul;
% dx = 0.5*mul;
dx = 0.05;

y0 = 2;
%% Continuation
Copt = struct('Nmax', 60, 'arclengthparm', 'orthogonal', ...
    'itDisplay', false, 'lsrch', 1, 'predictororder', 1, ...
    'DynScale', 1, 'Dscale', [1; 10], 'dsmin', 0.01, 'stepadapt', 0);
% Copt  =struct('Nmax', 50, 'dsmin', 0.01);
YXs = PRECOCONT(@(yx) PARABFUN(yx, c1, c2, po), y0, xstart, xend, dx, Copt);
% YXs = PRECOCONT_SC(@(yx) PARABFUN(yx, c1, c2, po), y0, xstart, xend, dx, Copt);
% YXs = CONTINUE(@(yx) PARABFUN(yx, c1, c2, po), y0, xstart, xend, dx, Copt);

%% plot
ys = sqrt(c1/c2)*linspace(-1, 1, 100);
switch po
    case 1
        xs = (c1-c2*ys.^2);
    case 2
        xs = (c1-c2*ys.^2).*ys;
end

figure(1)
clf()
plot(YXs(2,:), YXs(1,:), 'o-'); hold on
plot(xs, ys, 'k--')
grid on

% %% NLVIB
% Sopt = struct('Nmax', 1000, 'arclengthparm', 'normal', 'jac', 'off', 'dynamicDscale', 1);
% YXn = solve_and_continue(y0, @(ul) PARABFUN(ul, c1, c2), xstart, xend, dx, Sopt);
% 
% plot(YXn(end,:), YXn(1,:), 'x-')

%% 
function [R, dRdu, dRdl] = PARABFUN(ul, c1, c2, po)
    lam = ul(2);
    u = ul(1);
    
    switch po
        case 1 % Parabola
            R = lam-c1+c2*u^2;
            dRdu = 2*c2*u;
            dRdl = 1.0;
        case 2 % Cubic
            R = lam-(c1-c2*u^2)*u;
            dRdu = -c1+3*c2*u^2;
            dRdl = 1.0;
    end
%     dRdu = [dRdu dRdl];
end