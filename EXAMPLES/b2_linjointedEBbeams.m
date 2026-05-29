clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a linear jointed Euler-Bernoulli Beam example

%% Setup Model
Ey = 190e9;
rho = 7680;
wid = 0.2;  % Width
brd = 0.4;  % Breadth
Ar = wid*brd;  % Area
Iy = wid^3*brd/12;  % 2nd moment of area
L0 = 2.0;  % Total Length
Klib = struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25));
wcomps = [1 1;  % First component -> exp(  k x )
         -1 1;  % Second component-> exp( -k x )
         1j 1;  % Third component -> exp( ik x )
        -1j 1]; % Fourth component-> exp(-ik x )
% Setup "wave-based pieces"
pcs = [struct('coords', [0;L0/3;L0], 'wcomps', wcomps);
    struct('coords', [L0;2*L0], 'wcomps', wcomps)];
% Setup Boundary Conditions. Fix-Fix used here.
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j]);
    struct('i', 5, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j])];

%% Setup the linear Joint
kJ = 1e9;
cJ = 320;
cofs = @(w,xi) [(-Ey*Iy*Klib.K(w,xi)^3)*[1 -1 -1j 1j 0 0 0 0] +...
    (kJ-1j*cJ*w)*[1 1 1 1 -1 -1 -1 -1];
    (-Ey*Iy*Klib.K(w,xi)^2)*[1 1 -1 -1 0 0 0 0] +...
    (kJ-1j*cJ*w)*[1 -1 1j -1j -1 1 -1j 1j];
    (-Ey*Iy*Klib.K(w,xi)^3)*[1 -1 -1j 1j -1 1 1j -1j];
    (-Ey*Iy*Klib.K(w,xi)^2)*[1 1 -1 -1 -1 -1 1 1]];

joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs);

%% Setup Excitation
Mx = @(w,xi) inv([Ey*Iy*Klib.K(w,xi)^3*[-1 1 1j -1j];
                  Ey*Iy*Klib.K(w,xi)^2*[-1 -1 1 1];
                  [1 1 1 1];
                  Klib.K(w,xi)*[1 -1 1j -1j]]);

excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0], ...
    'rcofs0', [1/2;0;0;0]);
%'nh' sets the harmonic at which to apply the excitation

%% Pre-Processing
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);
Nwc = size(wcomps,1);  % Number of wave components

%% Conduct Linear Forced Response Analysis
Nw = 1000;
Ws = linspace(0, 2e3, Nw);
Npts = pcs(end).irange(end);  
ACs = zeros(Npts*Nwc,Nw);
Famp = 2e3;
for iw=1:Nw
    [Amat, ~, ~, Fv] = WVAMAT([Ws(iw);0], 1, pcs, bcs, joints, Klib);
    ACs(:,iw) = Amat\(Fv*Famp);
end

%% Plot Forced Response - EB
opi = 5:8; 
figure()
clf()

subplot(2,1,1)
semilogy(Ws/1e3, abs(2*sum(ACs(opi,:))));
grid on
ylabel('Response (m)')
hold on  

subplot(2,1,2)
plot(Ws/1e3, rad2deg(angle(2*sum(ACs(opi,:)))));
grid on
set(gca, 'YTick', -180:90:180)
ylabel('Phase (degs)')
xlabel('Frequency (k rad/s)')
hold on 

