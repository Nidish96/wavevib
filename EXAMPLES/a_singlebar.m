clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This shows the most basic example of a 1D bar fixed at both
%ends. 

%% Setup model
Ey = 2.62e11;
rho = 1280;
ell = 1.0;
Klib = struct('K', @(w,xi) w/sqrt(Ey/rho));
wcomps = [1j 1;  % First component -> exp( j k x )
         -1j 1]; % Sec component   -> exp(-j k x )

% Setup "wave-based pieces"
pcs = struct('coords', [0;ell], 'wcomps', wcomps);
% Note: The coordinates can be 3D also. The rows are interpreted as
% individual points. Coordinates not given are assumed to be zeros.

% Setup Boundary Conditions
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1]); % fixed end at point 1
    struct('i', 2, 'cofs', @(w,xi) [1 1])]; % fixed end at point 2
% Note: The BC coefficients are given as function handles in two
% parameters, (w,xi). The first is interpreted as frequency and the second
% is left to the discretion of the user. The number of columns has to be
% equal to the number of wave components in each case. 

%% Preprocess inputs
% This step computes the analytical derivatives for expressions given above
[pcs, bcs, ~, ~, Klib] = WBPREPROC(pcs, bcs, [], [], Klib);
% Note: the empty inputs are used to specify joints and excitation (see the
% other examples)

%% Compute determinant of linear Jacobian
Nw = 5000;
Ws = linspace(0, 1e6, Nw);
Ds = zeros(1,Nw);
for iw=1:Nw
    Ds(iw) = WVLDETFUN([Ws(iw);0], 1, pcs, bcs, [], Klib);
end

%% Plot
figure(1)
clf()
semilogy(Ws, Ds, '-')
xlabel('Frequncy (rad/s)')
ylabel('Jacobian Determinant')

%% Plot Highest Mode Shape
[pks, locs] = findpeaks(-Ds);
Wpks = Ws(locs);
Ks = Klib.K(Wpks,Wpks*0);

Nx = 100;
Xs = linspace(0, ell, Nx)';

mi = 20;
figure(2)
clf()
plot(Xs, sin(Ks(mi).*Xs), '.-');