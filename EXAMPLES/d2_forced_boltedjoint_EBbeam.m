clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a EB beam with single lap bolted joint under 
% forced vib. example 

%% Setup
Ey = 210e9;
rho = 7850;
wid = 1e-2;  % Width
brd = 1e-2;  % Breadth
Ar = wid*brd;  % Area
Iy = wid^3*brd/12;  % 2nd moment of area
L0 = 0.5;  % Total Length
alpha =0.25; %mass proportional damping coeff
beta = 2e-6; %stiffness proportional damping coeff

%Wavenumbers with mass & stiffness proportional damping terms
Klib = [struct('K', @(w,xi) ((rho*Ar/Ey/Iy) ...
    .*(w.^2 .*(1 - 1i*alpha./w)./(1 + 1i*beta.*w))).^(1/4) );
    struct('K', @(w,xi) sqrt((rho/Ey) ...
    .*(w.^2 .*(1 - 1i*alpha./w)./(1 + 1i*beta.*w))) )];

%Wave components
wcomps = [1j 1; % First component -> exp(ik1 x )
         -1j 1; % Second component-> exp(-ik1 x )
           1 1; % Third component -> exp(k1 x )
          -1 1; % Fourth component-> exp(-k1 x )
          1j 2; % Fifth component-> exp(ik2 x )
         -1j 2]; % Sixth component-> exp(-ik2 x )

pcs = [struct('coords',[0;L0/3;L0/2],'wcomps', wcomps);
    struct('coords',[L0/2;L0],'wcomps', wcomps)];  

%Fixed-Fixed BCs
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1 0 0; ...
    1j -1j 1 -1 0 0; 0 0 0 0 1 1]);
    struct('i', 5, 'cofs', @(w,xi) [1 1 1 1 0 0; ...
    1j -1j 1 -1 0 0; 0 0 0 0 1 1])];

%% Lap Jointed setup
lj = 0.03; %lap length
dia = 0.006 ; %diameter of the bolt
Eb = Ey ; %young's modulus of the bolt
n = 1 ; %no. of shear planes
a = 2/3; %Huth's Empricial constant
b = 3; %Huth's Empricial constant
h1 = brd;
h2 = brd;
E1 = Ey;
E2 = Ey;

%axial, transverse, and rotational stiffness 
k_x = 1/((((h1+h2)/(2*dia))^a) * (b/n) * ((1/(h1*E1)) + (1/(n*h2*E2)) ...
    + (1/(n*h1*Eb)) + (1/(2*n*h2*Eb)))); %Huth's Formula
k_bolt = (Eb * pi * (dia/2)^2) /(h1+h2);
k_comp = 0.5774*pi*Eb*dia/(2*log(5*((0.5774*(h1+h2)) ...
    +(0.5*dia))/((0.5774*(h1+h2))+(2.5*dia))));
k_z = k_comp + k_bolt;
k_theta = k_z*lj^2 / 12;
% damping coeff (cx & cz)
c_x = 3.645e5;
c_z = 3.645e5;

%Linear spring-damper Interaction forces and moments, Joint Equilibiriums 
cofs = @(w,xi)[Ey*Ar*Klib(2).K(w,xi)*[0 0 0 0 0 0 0 0 0 0 1j -1j] ...
    - k_x*[0 0 0 0 -1 -1 0 0 0 0 1 1] ...
    - (w*c_x)*[0 0 0 0 1j 1j 0 0 0 0 -1j -1j]; %F2-kx(u2-u1)-cx(u2dot-u1dot)=0
    (-Ey*Iy*Klib(1).K(w,xi)^3)*[0 0 0 0 0 0 -1j 1j 1 -1 0 0] ...
    - k_z*[-1 -1 -1 -1 0 0 1 1 1 1 0 0] ...
    - (w*c_z)*[1j 1j 1j 1j 0 0 -1j -1j -1j -1j 0 0]; %V2-kz(y2-y1)-cz(y2dot-y1dot)=0
    (-Ey*Iy*Klib(1).K(w,xi)^2)*[0 0 0 0 0 0 -1 -1 1 1 0 0] ...
    - k_theta*[-1j 1j -1 1 0 0 1j -1j 1 -1 0 0]; %M2-ktheta(psi2-psi1)=0
    (-Ey*Iy*Klib(1).K(w,xi)^3)*[1j -1j -1 1 0 0 -1j 1j 1 -1 0 0]; %V2-V1=0
    [0 0 0 0 -1j 1j 0 0 0 0 1j -1j]; %F2-F1=0   
    (-Ey*Iy*Klib(1).K(w,xi)^2)*[ 1 1 -1 -1 0 0 -1 -1 1 1 0 0] ...
    + (-0.5*h1*Ey*Ar*Klib(2).K(w,xi))*[0 0 0 0 1j -1j 0 0 0 0 1j -1j] + ...
    (-0.5*dia*Ey*Iy*Klib(1).K(w,xi)^3)*[-1j 1j 1 -1 0 0 -1j 1j 1 -1 0 0] 
    %M2-M1-h/2(F1+F2)+d/2(V1+V2)=0
    ];

joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs);

%% Excitation
Mx = @(w,xi)inv([(-Ey*Iy*Klib(1).K(w,xi)^3)*[-1j 1j 1 -1 0 0];
                 (-Ey*Iy*Klib(1).K(w,xi)^2)*[-1 -1 1 1 0 0];
                 (Ey*Ar*Klib(2).K(w,xi))*[0 0 0 0 1j -1j];
                  [1 1 1 1 0 0];
                  [1j -1j 1 -1 0 0];
                  [0 0 0 0 1 1]]);

excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0;0;0], ...
    'rcofs0', [1/2;0;0;0;0;0]);

%% Preprocess Everything
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);
Nwc = size(wcomps,1);  % Number of wave components

%% Conduct Linear Forced Response Analysis
Nw = 1000;
Ws = linspace(0, 1.5e4, Nw);
Npts = pcs(end).irange(end);  
ACs = zeros(Npts*Nwc,Nw);
Famp = 1e2;
for iw=1:Nw
    [Amat, ~, ~, Fv] = WVAMAT([Ws(iw);0], 1, pcs, bcs, joints, Klib);
    ACs(:,iw) = Amat\(Fv*Famp);
end

%% Plot Forced Response
opi = 7:10;  % Output wave coefficients
figure()
clf;
subplot(2,1,1)
semilogy(Ws/1e3, abs(2*sum(ACs(opi,:))));
grid on
ylabel('Response (m)')
xlabel('Frequency (k rad/s)')
subplot(2,1,2)
plot(Ws/1e3, rad2deg((angle(2*sum(ACs(opi,:))))));
grid on
set(gca, 'YTick', -180:90:180)
ylabel('Phase (degs)')
xlabel('Frequency (k rad/s)')

