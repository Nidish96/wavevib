clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This shows the linear forced response computation for the 1D
%bar example considered in Balaji, Brake, Leamy (2022a,2022b)

%% Setup model
Ey = 2.62e11;
rho = 1280;
ell = 1.0;
Ar = 31.75e-3*54e-3;
Klib = struct('K', @(w,xi) w/sqrt(Ey/rho));
wcomps = [1j 1;  % First component -> exp( j k x )
         -1j 1]; % Sec component   -> exp(-j k x )

% Setup "wave-based pieces"
pcs = [struct('coords', [0;0.28;ell/3], 'wcomps', wcomps);
    struct('coords', [ell/3;ell], 'wcomps', wcomps)];

% Setup Boundary Conditions (pts indexed as given above in sequence)
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1]);
    struct('i', 5, 'cofs', @(w,xi) [1 1])];

%% Setup the linear Joint
kJ = 1e9;
cJ = 320;
cofs = @(w,xi) [(1j*Klib.K(w,xi)*Ey*Ar)*[1, -1, 0, 0] +...
    (kJ-1j*cJ*w)*[1 1 -1 -1];1, -1, -1, 1];
% Coefficients are to be supplied such that 
% COFS*a = 0    represents the joint, where
% a is the vector of the wave coefficients of the two pts.
joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs);
%NOTE: The 'type' parameter specifies how many points are joined at this
%location. So far only type 2 has been implemented.

%% Setup Excitation
excs = struct('i', 2, ...
    'rcofs', @(w,xi) (1/2/(2j*Klib.K(w,xi)*Ey*Ar))*[-1;1]);
%'i' specifies the point of excitation. 
%'rcofs' specifies the coefficients in the RHS of the eqn such that,
%   [I -I]*[a;b] = rcofs        represents the excitation.
%   where "a" and "b" are the vector of wave coefficients just before and
%   after the point of excitation.

%% Preprocess Everything
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);
Nwc = size(wcomps,1);  % Number of wave components

%% Conduct Linear Forced Response Analysis
Nw = 1000;
Ws = linspace(0, 1e5, Nw);
Npts = pcs(end).irange(end);  % Total number of points. 
                              % Note: This is NOT the same as before. 
                              % Preprocessing usually adds more/reorders points.
ACs = zeros(Npts*Nwc,Nw);

Famp = 150e5;  % 15MN Excitation amplitude
for iw=1:Nw
    % The following function calculates the "A" matrix and the "Force"
    % vector (both complex). An argument 'r' can be passed as the last
    % argument to obtain everything in real form (real and imaginary parts 
    % stored separately). "Fv" is the force vector assuming unit
    % excitation.
    [Amat, ~, ~, Fv] = WVAMAT([Ws(iw);0], 1, pcs, bcs, joints, Klib);
    ACs(:,iw) = Amat\(Fv*Famp);
end

%% Plot Forced Response
opi = 9:10;  % Output wave coefficients
figure(1)
clf()
subplot(2,1,1)
semilogy(Ws/1e3, abs(2*sum(ACs(opi,:))));
grid on
ylabel('Response (m)')
subplot(2,1,2)
plot(Ws/1e3, rad2deg(angle(2*sum(ACs(opi,:)))));
grid on
set(gca, 'YTick', -180:90:180)
ylabel('Phase (degs)')
xlabel('Frequency (k rad/s)')