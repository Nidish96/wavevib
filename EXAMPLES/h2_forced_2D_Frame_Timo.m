clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is forced vibration of a portal frame with Timoshenko members example
% Initially, the beam-column joints are treated as rigid angle joints.
% In Nonlinear example, they are replaced with a nonlinear spring damper joint.
% Parameters are taken from M. I. Basci et al, 1978.

%% Setup Model
Ey = 2.068e11;
nu = 0.3;
G = Ey/(2*(1+nu)); %shear modulus
kappa = 10*(1+nu)/(12+(11*nu)); %shear correction factor
rho = 7842.22747;
wid = 0.89208;
brd = 0.021255;
Ar = 0.0189612524;  % Area
Iy = 0.001257435137946;  % 2nd moment of area
% Coefficients
Cb = sqrt(Ey*Iy / (rho*Ar)); %bending stiffness
Cs = sqrt(kappa*G*Ar/(rho*Ar)); %shear stiffness
Cr = sqrt(Iy / Ar); %rotational effects
% Klib K1, K2 & K3
 Klib = [struct('K', @(w,xi) sqrt( ...
     0.5*((1/Cs)^2 + (Cr/Cb)^2)*w.^2 + ...
     sqrt((w.^2)/(Cb^2) + 0.25*((1/Cs)^2 - (Cr/Cb)^2)^2 * w.^4) ));
     struct('K', @(w,xi) sqrt( ...
     abs(0.5*((1/Cs)^2 + (Cr/Cb)^2)*w.^2 - ...
     sqrt((w.^2)/(Cb^2) + 0.25*((1/Cs)^2 - (Cr/Cb)^2)^2 * w.^4))));
     struct('K', @(w,xi) w*sqrt(rho/Ey))];
% Wave components
wcomps = [-1j 1; % First component -> exp(-ik1 x )
           -1 2; % Second component-> exp(-k2 x )
           1j 1; % Third component -> exp(ik1 x )
            1 2; % Fourth component-> exp(k2 x )
          -1j 3; % Fifth component-> exp(-ik3 x )
          1j 3]; % Sixth component-> exp(ik3 x )
% The relations between the coefficients of deflection wave components 
% & bending slope wave components
P = @(w,xi) (Klib(1).K(w,xi).^2 * Cs^2 - w.^2) ./ (Klib(1).K(w,xi) * Cs^2);
N = @(w,xi) (Klib(2).K(w,xi).^2 * Cs^2 + w.^2) ./ (Klib(2).K(w,xi) * Cs^2);
% Pieces
L = 15.24;  
H = 15.24;
pcs = [struct('coords',[0 0;0 H/2;0 H],'wcomps', wcomps);
    struct('coords',[0 H;L H],'wcomps', wcomps);   
    struct('coords',[L H;L 0],'wcomps', wcomps)]; 
% BCs (Fix_Fix)
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1 0 0; ...
    -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0; 0 0 0 0 1 1]);
    struct('i', 7, 'cofs', @(w,xi) [1 1 1 1 0 0; ...
    -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0; 0 0 0 0 1 1])];

%% Angle Joints
m = brd*brd*wid*rho;
J = (brd^2 + brd^2)*m/12;

cofs = @(w,xi)[ [0 0 0 0 1 1 0 0 0 0 0 0] - [0 0 0 0 0 0 1 1 1 1 0 0];%u1-v2=0
    [0 0 0 0 0 0 0 0 0 0 1 1] + [1 1 1 1 0 0 0 0 0 0 0 0];%u2+v1=0
    [-1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0 ...
    1j*P(w,xi) 1*N(w,xi) -1j*P(w,xi) -1*N(w,xi) 0 0];% psi1-psi2=0
    (G*Ar*kappa)*[0 0 0 0 0 0 -1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) ...
    1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0] ...
    - (Ey*Ar*Klib(3).K(w,xi))*[0 0 0 0 -1j 1j 0 0 0 0 0 0] ...
    + (m*w*w)*[0 0 0 0 1 1 0 0 0 0 0 0]; % V2-F1=0
    (G*Ar*kappa)*[-1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) ...
    1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0 0 0 0 0 0 0] ...
    + (Ey*Ar*Klib(3).K(w,xi))*[0 0 0 0 0 0 0 0 0 0 -1j 1j] ...
    + (m*w*w)*[0 0 0 0 0 0 0 0 0 0 1 1]; % F2+V1=0
    (Ey*Iy)*[(Klib(1).K(w,xi)*P(w,xi)) (-Klib(2).K(w,xi)*N(w,xi)) (Klib(1).K(w,xi)*P(w,xi)) ...
    (-Klib(2).K(w,xi)*N(w,xi)) 0 0 (-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) ...
    (-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) 0 0] ...
    + (0.5*brd*G*Ar*kappa)*[-1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) ...
    1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0 -1j*(Klib(1).K(w,xi)-P(w,xi)) ...
    (-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0] ...
    + (J*w*w)*[-1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0 0 0 0 0 0 0]]; 
    %M2-M1+(h/2 * (V1+V2)) = J(w^2)psij
  
joints = [struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs);
          struct('type', 2, 'i', 5, 'j', 6, 'cofs', cofs)];

%% Excitation
Mx = @(w,xi)inv([G*Ar*kappa*[-1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) ...
    1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0];
    (Ey*Iy)*[(-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) ...
    (-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) 0 0];
    Ey*Ar*Klib(3).K(w,xi)*[0 0 0 0 -1j 1j];
    [1 1 1 1 0 0];
    [-1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0];
    [0 0 0 0 1 1]]);

excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0;0;0], ...
    'rcofs0', [1/2;0;0;0;0;0]);

%% Preprocess Everything
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);
Nwc = size(wcomps,1);  % Number of wave components

%% Conduct Linear Forced Response Analysis
Nw = 500;
Ws = linspace(1, 2e2, Nw);
Npts = pcs(end).irange(end);  
ACs = zeros(Npts*Nwc,Nw);
Famp = 100;
for iw=1:Nw
    [Amat, ~, ~, Fv] = WVAMAT([Ws(iw);0], 1, pcs, bcs, joints, Klib);
    ACs(:,iw) = Amat\(Fv*Famp);
end

%% Plot Forced Response
opi = 7:12;  % Output wave coefficients
figure()
clf;
subplot(2,1,1)
semilogy(Ws/1e3, abs(2*sum(ACs(opi,:))));
grid on
ylabel('Response (m)')
subplot(2,1,2)
plot(Ws/1e3, rad2deg(unwrap(angle(2*sum(ACs(opi,:))))));
grid on
%set(gca, 'YTick', -180:90:180)
ylabel('Phase (degs)')
xlabel('Frequency (k rad/s)')