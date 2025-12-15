clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is free vibration of a portal frame with EB members example 
% Parameters are taken from M. I. Basci et al, 1978.
%% Setup Model
Ey = 2.068e11;
rho = 7842.22747;
wid = 0.89208;
brd = 0.02125;
Ar = 0.0189612524;  % Area
Iy = 0.001257435137946;  % 2nd moment of area

Klib = [struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25));
        struct('K', @(w,xi) w*sqrt(rho/Ey))];
%Wave components
wcomps = [1 1;  % -> exp(k1 x )
         -1 1;  % -> exp(-k1 x )
         1j 1;  % -> exp(ik1 x )
        -1j 1;  % -> exp(-ik1 x )
         1j 2;  % -> exp(ik2 x )
        -1j 2]; % -> exp(-ik2 x )
%pieces
L = 15.24;  
H = 15.24; 
pcs = [struct('coords',[0 0;0 H],'wcomps', wcomps);
       struct('coords',[0 H;L H],'wcomps', wcomps);   
       struct('coords',[L H;L 0],'wcomps', wcomps)]; 
% BCs (Fix_Fix)
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1 0 0; 1 -1 1j -1j 0 0; 0 0 0 0 1 1]);
       struct('i', 6, 'cofs', @(w,xi) [1 1 1 1 0 0; 1 -1 1j -1j 0 0; 0 0 0 0 1 1])];


%joint
m =  brd*brd*wid*rho;
J = (brd^2 + brd^2)*m/12;
cofs3 = @(w,xi)[ [0 0 0 0 1 1 0 0 0 0 0 0] - [0 0 0 0 0 0 1 1 1 1 0 0];%u1-v2=0
    [0 0 0 0 0 0 0 0 0 0 1 1] + [1 1 1 1 0 0 0 0 0 0 0 0];%u2+v1=0
    [1 -1 1j -1j 0 0 -1 1 -1j 1j 0 0];% psi1-psi2=0
    (-Ey*Iy*Klib(1).K(w,xi)^3)*[0 0 0 0 0 0 1 -1 -1j 1j 0 0] - (Ey*Ar*Klib(2).K(w,xi))*[0 0 0 0 1j -1j 0 0 0 0 0 0] + (m*w.^2)*[0 0 0 0 1 1 0 0 0 0 0 0]; % V2-F1=0
    (Ey*Ar*Klib(2).K(w,xi))*[0 0 0 0 0 0 0 0 0 0 1j -1j] + (-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0 0 0 0 0 0 0] + (m*w.^2)*[0 0 0 0 0 0 0 0 0 0 1 1]; % F2+V1=0
    (Ey*Iy*Klib(1).K(w,xi)^2)*[-1 -1 1 1 0 0 1 1 -1 -1 0 0] + (-0.5*brd*Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0 1 -1 -1j 1j 0 0] + (J*w.^2)*[1 -1 1j -1j 0 0 0 0 0 0 0 0]]; % M2-M1+(h/2 * (V1+V2)) = 0]

joints = [struct('type', 2, 'i', 2, 'j', 3, 'cofs', cofs3);
          struct('type', 2, 'i', 4, 'j', 5, 'cofs', cofs3)];

%% Pre-Processing
[pcs, bcs, joints, ~, Klib] = WBPREPROC(pcs, bcs, joints, [], Klib);

%% 
Nw = 1000;
Ws = linspace(1, 1.4e2, Nw);
Ds = zeros(1,Nw);
for iw=1:Nw
    Ds(iw) = WVLDETFUN([Ws(iw);0], 1, pcs, bcs, joints, Klib);
end
% WVLDETFUN([Ws(iw);0], 1, pcs, bcs, joints, Klib);

%% Plot
figure(1)
clf()
semilogy(Ws, Ds, '-')
xlabel('Frequncy (rad/s)')
ylabel('Jacobian Determinant')
