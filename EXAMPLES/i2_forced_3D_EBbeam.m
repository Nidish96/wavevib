clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a forced vibration of single 3D/spatial Euler-Bernoulli Beam example
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
Klib = [struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25)); 
        %wavenumber for bending in y direction
        struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iz)^(0.25)); 
        %wavenumber for bending in z direction
        struct('K', @(w,xi) w*sqrt(rho/Ey))]; %longitudinal wavenumber

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

% Setup "wave-based pieces"
pcs = struct('coords', [0 0 0;L0/3 0 0;L0 0 0], 'wcomps', wcomps);
% Setup Boundary Conditions. Fix-Fix used here.
bcs = [struct('i', 1, 'cofs', @(w,xi) [0 0 0 0 0 0 0 0 1 1; ...
    1 1 1 1 0 0 0 0 0 0; 0 0 0 0 1 1 1 1 0 0; ...
    1 -1 1j -1j 0 0 0 0 0 0;0 0 0 0 1 -1 1j -1j 0 0]);
    struct('i', 3, 'cofs', @(w,xi) [0 0 0 0 0 0 0 0 1 1; ...
    1 1 1 1 0 0 0 0 0 0; 0 0 0 0 1 1 1 1 0 0;1 -1 1j -1j 0 0 0 0 0 0; ...
    0 0 0 0 1 -1 1j -1j 0 0])];

%% Excitation
Mx = @(w,xi)inv([(-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0 0 0 0 0];
    (-Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0 0 0 0 0];
    [1 1 1 1 0 0 0 0 0 0];
    Klib(1).K(w,xi)*[1 -1 1j -1j 0 0 0 0 0 0];
    (-Ey*Iz*Klib(2).K(w,xi)^3)*[0 0 0 0 1 -1 -1j 1j 0 0];
    (-Ey*Iz*Klib(2).K(w,xi)^2)*[0 0 0 0 1 1 -1 -1 0 0];
    [0 0 0 0 1 1 1 1 0 0];
    Klib(2).K(w,xi)*[0 0 0 0 1 -1 1j -1j 0 0];
    (Ey*Ar*Klib(3).K(w,xi))*[0 0 0 0 0 0 0 0 1j -1j];
    [0 0 0 0 0 0 0 0 1 1]]);
excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0;1/2;0;0;0;0;0], ...
    'rcofs0', [1/2;0;0;0;1/2;0;0;0;0;0]);  %Excitation in y and z direction to excite all modes

%% Setup the linear Joint
kJ = 1e9;
cJ = 320;
cofs = @(w,xi) [(-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0] +...
    (kJ-1j*cJ*w)*[1 1 1 1 0 0 0 0 0 0 -1 -1 -1 -1 0 0 0 0 0 0];
    (-Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0] +...
    (kJ-1j*cJ*w)*[1 -1 1j -1j 0 0 0 0 0 0 -1 1 -1j 1j 0 0 0 0 0 0];
    (-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j -1 1 1j -1j];
    (-Ey*Iy*Klib(2).K(w,xi)^2)*[1 1 -1 -1 -1 -1 1 1];
    (-Ey*Iz*Klib(2).K(w,xi)^3)*[0 0 0 0 0 0 1 -1 -1j 1j 0 0 0 0 0 0 0 0 0 0] +...
    (kJ-1j*cJ*w)*[0 0 0 0 0 0 1 1 1 1 0 0 0 0 -1 -1 -1 -1 0 0 0 0 0 0];
    (-Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0] +...
    (kJ-1j*cJ*w)*[1 -1 1j -1j 0 0 0 0 0 0 -1 1 -1j 1j 0 0 0 0 0 0];
    (-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j -1 1 1j -1j];
    (-Ey*Iy*Klib(2).K(w,xi)^2)*[1 1 -1 -1 -1 -1 1 1]
    ];

joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs);

%% Pre-Processing
[pcs, bcs, ~, excs, Klib] = WBPREPROC(pcs, bcs, [], excs, Klib);
Nwc = size(wcomps,1);  % Number of wave components

%% Conduct Linear Forced Response Analysis
Nw = 1000;
Ws = linspace(0, 8e2, Nw);
Npts = pcs(end).irange(end);  
ACs = zeros(Npts*Nwc,Nw);
Famp = 1e2;
for iw=1:Nw
    [Amat, ~, ~, Fv] = WVAMAT([Ws(iw);0], 1, pcs, bcs, [], Klib);
    ACs(:,iw) = Amat\(Fv*Famp);
end

%% Plot Forced Response
opi = 11:20; 
figure()
clf()
subplot(2,1,1)
semilogy(Ws/1e3, abs(2*sum(ACs(opi,:))));
grid on
ylabel('Response (m)')
hold on  
subplot(2,1,2)
plot(Ws/1e3, rad2deg(unwrap(angle(2*sum(ACs(opi,:))))));
grid on
set(gca, 'YTick', -180:90:180)
ylabel('Phase (degs)')
xlabel('Frequency (k rad/s)')
hold on 



