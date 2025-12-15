clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

addpath('../ROUTINES/export_fig')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a nonlinear jointed Euler-Bernoulli Beam example
%similar to the jointed bar example. Parameters taken from 
% Krishna and Chandramouli, 2012

savfig = false;
animfig = false;
savdat = true;
analyze = true;
%% Setup Model
Ey = 2.1e11;
rho = 7680;
thk = 3e-3;   % Thickness 1
wid  = 40e-3;  % Width
Ar = thk*wid;  % Area
Iy = thk^3*wid/12;  % 2nd moment of area
L1  = 560e-3;  % Primary beam length
xF1 = L1/6;    % Excitation Location (on primary beam)
L2 = 445e-3;  % Secondary beam length (each)
% Klib = struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25));  % Undamped EB-Beam Case
al = 0.80;    % 0.80, 1.80
bt = 1.1e-4;  % 1.1e-4, 2.475e-4
% Joint parameters
knl = 1210/2;      % 220,      880,    (242, 484, 880, 1210)
gap = 2.5e-3;   % 2.5e-3,   0.35,   0.35

Klib = struct('K', @(w,xi) ((rho*Ar*(w.^2+1j*w*al))./(Ey*Iy*(1-1j*w*bt))).^(0.25) );
wcomps = [1 1;  % First component -> exp(  k x )
         -1 1;  % Second component-> exp( -k x )
         1j 1;  % Third component -> exp( ik x )
        -1j 1]; % Fourth component-> exp(-ik x )

pcs = struct('coords', [0; 0.5; 1.0], 'wcomps', wcomps);
bcs = [struct('i', 1, 'cofs', @(w, xi) [1 1 1 1;1 1 -1 -1], ...
              'cofs0', @(xi) [1 0 0 0;0 0 1 0]);
       struct('i', 3, 'cofs', @(w, xi) [1 1 1 1;1 1 -1 -1], ...
              'cofs0', @(xi) [1 0 0 0;0 0 1 0])];
Mx = @(w,xi) inv([Ey*Iy*Klib.K(w,xi)^3*[-1 1 1j -1j];
                  Ey*Iy*Klib.K(w,xi)^2*[-1 -1 1 1];
                  [1 1 1 1];
                  Klib.K(w,xi)*[1 -1 1j -1j]]);
excs = struct('i', 2, 'nh', 0, ...
              'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0], ...
              'rcofs0', @(xi) [0;0;0;1/(6*Ey*Iy)]);

[pcs, bcs, ~, excs, Klib] = WBPREPROC(pcs, bcs, [], excs, Klib);

%%
Wst = 0.5;

h = 0;
[Amat, ~, ~, Fv] = WVAMAT([Wst;0], h, pcs, bcs, [], Klib);
