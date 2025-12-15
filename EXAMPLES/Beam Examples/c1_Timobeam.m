clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a single Timoshenko Beam example similar to Euler
%Bernoulli Beam example
%% Setup Model
Ey = 190e9;
nu = 0.29;
G = 73.6e9; %shear modulus
kappa = 10*(1+nu)/(12+(11*nu)); %shear correction factor
rho = 7680;
wid = 0.4;  % Width
brd = 0.4;  % Breadth
Ar = wid*brd;  % Area
Iy = wid^3*brd/12;  % 2nd moment of area
L0 = 4.0;  % Total Length

%coefficients
Cb = sqrt(Ey*Iy / (rho*Ar)); %bending stiffness
Cs = sqrt(kappa*G*Ar/(rho*Ar)); %shear stiffness+
Cr = sqrt(Iy / Ar); %rotational effects

%Wavenumbers K1 & K2
Klib = [struct('K', @(w,xi) sqrt( ...
     0.5*((1/Cs)^2 + (Cr/Cb)^2)*w.^2 + ...
     sqrt((w.^2)/(Cb^2) + 0.25*((1/Cs)^2 - (Cr/Cb)^2)^2 * w.^4) ));
     struct('K', @(w,xi) sqrt( ...
     (-0.5*((1/Cs)^2 + (Cr/Cb)^2)*w.^2 + ...
     sqrt((w.^2)/(Cb^2) + 0.25*((1/Cs)^2 - (Cr/Cb)^2)^2 * w.^4))))];

%Wave components
wcomps = [1j 1;  % First component -> exp(ik1 x )
         1 2;  % Second component-> exp(k2 x )
         -1j 1;  % Third component -> exp(-ik1 x )
        -1 2]; % Fourth component-> exp(-k2 x )

%pieces
pcs = [struct('coords', [0;L0], 'wcomps', wcomps)];

% The relations between the coefficients of deflection wave components & bending slope wave components
P = @(w,xi) (Klib(1).K(w,xi))*( 1  - ((w.^2) ./ ((Klib(1).K(w,xi).^2) * Cs^2)));
N = @(w,xi) (Klib(2).K(w,xi))*( 1  + ((w.^2) ./ ((Klib(2).K(w,xi).^2) * Cs^2)));

% Fixed-Fixed Boundary conditions 
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1; 1j*P(w,xi) 1*N(w,xi) -1j*P(w,xi) -1*N(w,xi)]);
    struct('i', 2, 'cofs', @(w,xi) [1 1 1 1; 1j*P(w,xi) 1*N(w,xi) -1j*P(w,xi) -1*N(w,xi)])];

%% Pre-processing
[pcs, bcs, ~, ~, Klib] = WBPREPROC(pcs, bcs, [], [], Klib);

%% Compute determinant of linear Jacobian
Nw = 1000;
Ws = linspace(1, 1e3, Nw);
Ds = zeros(1,Nw);
for iw=1:Nw
    Ds(iw) = WVLDETFUN([Ws(iw);0], 1, pcs, bcs, [], Klib);
end

%% Plot
figure(1); clf();
semilogy(Ws, Ds);
hold on;
xlabel('Frequency (rad/s)')
ylabel('Jacobian Determinant')




