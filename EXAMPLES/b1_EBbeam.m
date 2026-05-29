clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a single Euler-Bernoulli Beam example

%% Setup Model
Ey = 190e9;
rho = 7680;
wid = 0.2;  % Width
brd = 0.4;  % Breadth
Ar = wid*brd;  % Area
Iy = wid^3*brd/12;  % 2nd moment of area
L0 = 4.0;  % Total Length
Klib = struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25));
wcomps = [1 1;  % First component -> exp(  k x )
         -1 1;  % Second component-> exp( -k x )
         1j 1;  % Third component -> exp( ik x )
        -1j 1]; % Fourth component-> exp(-ik x )
% Setup "wave-based pieces"
pcs = struct('coords', [0;L0], 'wcomps', wcomps);
% Setup Boundary Conditions. Fix-Fix used here.
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j]);
    struct('i', 2, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j])];

%% Pre-processing
[pcs, bcs, ~, ~, Klib] = WBPREPROC(pcs, bcs, [], [], Klib);

%% Compute determinant of linear Jacobian
Nw = 1000;
Ws = linspace(eps, 5e3, Nw);
Ds = zeros(1,Nw);
for iw=1:Nw
    Ds(iw) = WVLDETFUN([Ws(iw);0], 1, pcs, bcs, [], Klib);
end
%% Plot
figure(); clf();
semilogy(Ws, Ds);
hold on;
xlabel('Frequency (rad/s)')
ylabel('Jacobian Determinant')
