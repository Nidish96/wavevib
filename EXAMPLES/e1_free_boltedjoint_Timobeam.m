clc
clear all
addpath('D:\MTech_IITM_Project\wavevib\ROUTINES\SOLVERS')
addpath('D:\MTech_IITM_Project\wavevib\ROUTINES\WBM')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a Timoshenko beam with single lap bolted joint under 
% free vib. example

%% Setup
Ey = 210e9;
nu = 0.3;
G = 80e9; %shear modulus
kappa = 10*(1+nu)/(12+(11*nu)); %shear correction factor
rho = 7850;
wid = 1e-2;  % Width
brd = 1e-2;  % Breadth
Ar = wid*brd;  % Area
Iy = wid^3*brd/12;  % 2nd moment of area
L0 = 0.5;  % Total Length

%coefficients
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

%Wave components
wcomps = [1j 1; % First component -> exp(ik1 x )
         -1j 1; % Second component-> exp(-ik1 x )
           1 2; % Third component -> exp(k2 x )
          -1 2; % Fourth component-> exp(-k2 x )
          1j 3; % Fifth component-> exp(ik3 x )
         -1j 3]; % Sixth component-> exp(-ik3 x )

% The relations between the coefficients of deflection wave components
%  & bending slope wave components
P = @(w,xi) (Klib(1).K(w,xi).^2 * Cs^2 - w.^2) ./ (Klib(1).K(w,xi) * Cs^2);
N = @(w,xi) (Klib(2).K(w,xi).^2 * Cs^2 + w.^2) ./ (Klib(2).K(w,xi) * Cs^2);

pcs = [struct('coords',[0;L0/2],'wcomps', wcomps);
    struct('coords',[L0/2;L0],'wcomps', wcomps)];  

bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1 0 0; ...
    1j*P(w,xi) -1j*P(w,xi) 1*N(w,xi) -1*N(w,xi) 0 0; 0 0 0 0 1 1]);
    struct('i', 4, 'cofs', @(w,xi) [1 1 1 1 0 0; ...
    1j*P(w,xi) -1j*P(w,xi) 1*N(w,xi) -1*N(w,xi) 0 0; 0 0 0 0 1 1])];

%% Lap Jointed setup
lj = 0.03; %lap length
dia = 0.006 ; %diameter of the bolt
Eb= Ey ; %young's modulus of the bolt
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
cofs = @(w,xi)[Ey*Ar*Klib(3).K(w,xi)*[0 0 0 0 0 0 0 0 0 0 1j -1j] ...
    - k_x*[0 0 0 0 -1 -1 0 0 0 0 1 1] ...
    - (w*c_x)*[0 0 0 0 1j 1j 0 0 0 0 -1j -1j]; %F2-kx(u2-u1)-cx(u2dot-u1dot)=0
    G*Ar*kappa*[0 0 0 0 0 0 1j*(Klib(1).K(w,xi)-P(w,xi)) 1j*(-Klib(1).K(w,xi)+P(w,xi)) ...
    (Klib(2).K(w,xi)-N(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 0 0] ...
    - k_z*[-1 -1 -1 -1 0 0 1 1 1 1 0 0] ...
    - (w*c_z)*[1j 1j 1j 1j 0 0 -1j -1j -1j -1j 0 0]; %V2-kz(y2-y1)-cz(y2dot-y1dot)=0
    Ey*Iy*[0 0 0 0 0 0 (-Klib(1).K(w,xi)*P(w,xi)) (-Klib(1).K(w,xi)*P(w,xi)) ...
    (Klib(2).K(w,xi)*N(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) 0 0] ...
    - k_theta*[-1j*P(w,xi) 1j*P(w,xi) -1*N(w,xi) 1*N(w,xi) 0 0 1j*P(w,xi) ...
    -1j*P(w,xi) 1*N(w,xi) -1*N(w,xi) 0 0]; %M2-ktheta(psi2-psi1)=0
    [1j*(-Klib(1).K(w,xi)+P(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) ...
    (Klib(2).K(w,xi)-N(w,xi)) 0 0 1j*(Klib(1).K(w,xi)-P(w,xi)) 1j*(-Klib(1).K(w,xi)+P(w,xi)) ...
    (Klib(2).K(w,xi)-N(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 0 0]; %V2-V1=0
    [0 0 0 0 -1j 1j 0 0 0 0 1j -1j]; %F2-F1=0
    Ey*Iy*[(Klib(1).K(w,xi)*P(w,xi)) (Klib(1).K(w,xi)*P(w,xi)) (-Klib(2).K(w,xi)*N(w,xi)) ...
    (-Klib(2).K(w,xi)*N(w,xi)) 0 0 (-Klib(1).K(w,xi)*P(w,xi)) (-Klib(1).K(w,xi)*P(w,xi)) ...
    (Klib(2).K(w,xi)*N(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) 0 0] + ...
    (-0.5*h1*Ey*Ar*Klib(3).K(w,xi))*[0 0 0 0 1j -1j 0 0 0 0 1j -1j] + ...
    (0.5*dia*G*Ar*kappa)*[ 1j*(Klib(1).K(w,xi)-P(w,xi)) 1j*(-Klib(1).K(w,xi)+P(w,xi)) ...
    (Klib(2).K(w,xi)-N(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 0 0 1j*(Klib(1).K(w,xi)-P(w,xi)) ...
    1j*(-Klib(1).K(w,xi)+P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 0 0]; 
    %M2-M1-h/2(F1+F2)+d/2(V1+V2)=0
    ];

joints = struct('type', 2, 'i', 2, 'j', 3, 'cofs', cofs);

%% Pre-Processing
[pcs, bcs, joints, ~, Klib] = WBPREPROC(pcs, bcs, joints, [], Klib);

%% Compute determinant of linear Jacobian
Nw = 1000;
Ws = linspace(0, 1.5e4, Nw);
Ds = zeros(1,Nw);
for iw=1:Nw
    Ds(iw) = WVLDETFUN([Ws(iw);0], 1, pcs, bcs, joints, Klib);
end

%% Plot
figure()
clf()
semilogy(Ws, Ds, '-')
xlabel('Frequncy (rad/s)')
ylabel('Jacobian Determinant')