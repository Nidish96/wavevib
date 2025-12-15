clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a free vibration of single 3D/spatial Euler-Bernoulli Beam example
%Parameters are taken from Ljiljana, et al, 2016.
%% Setup Model
Ey = 3e10;
nu = 0.25;
G = Ey/(2*(1+nu));
rho = 2548.42;
wid = 0.3; 
brd = 0.4;  
Ar = wid*brd;  
Iy = wid^3*brd/12; 
Iz = brd^3*wid/12;  
L0 = 8.0; 
Klib = [struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25)); %wavenumber for bending in y direction
        struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iz)^(0.25)); %wavenumber for bending in z direction
        struct('K', @(w,xi) w*sqrt(rho/Ey))]; %torisonal wavenumber
wcomps = [1 1;  % exp(k1 x )
         -1 1;  % exp(-k1 x )
         1j 1;  % exp(ik x )
        -1j 1;  % exp(ik1 x)
          1 2;  % exp(k2 x )
         -1 2;  % exp(-k2 x )
         1j 2;  % exp(ik2 x )
        -1j 2;  % exp(-ik2 x )
         1j 3;  % exp(ik3 x )
        -1j 3]; % exp(-ik3 x )
% pieces
pcs = struct('coords', [0 0 0;L0 0 0], 'wcomps', wcomps);
% Boundary Conditions. Fix-Fix used here.
bcs = [struct('i', 1, 'cofs', @(w,xi) [0 0 0 0 0 0 0 0 1 1; 1 1 1 1 0 0 0 0 0 0; 0 0 0 0 1 1 1 1 0 0;1 -1 1j -1j 0 0 0 0 0 0 ;0 0 0 0 1 -1 1j -1j 0 0]);
    struct('i', 2, 'cofs', @(w,xi) [0 0 0 0 0 0 0 0 1 1; 1 1 1 1 0 0 0 0 0 0; 0 0 0 0 1 1 1 1 0 0;1 -1 1j -1j 0 0 0 0 0 0 ;0 0 0 0 1 -1 1j -1j 0 0])];
%% Pre-processing
[pcs, bcs, ~, ~, Klib] = WBPREPROC(pcs, bcs, [], [], Klib);
%% Compute determinant of linear Jacobian
Nw = 1000;
Ws = linspace(0, 8e2, Nw);
Ds = zeros(1,Nw);
for iw=1:Nw
    Ds(iw) = WVLDETFUN([Ws(iw);0], 1, pcs, bcs, [], Klib);
end
% Plot
figure(2); clf();
semilogy(Ws, Ds);
hold on;
xlabel('Frequency (rad/s)')
ylabel('Jacobian Determinant')

