clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%% Setup Model
Ey = 190e9;
nu = 0.29;
G = Ey/(2*(1+nu)); %shear modulus
kappa = 10*(1+nu)/(12+(11*nu)); %shear correction factor
rho = 7680;
wid = 0.2;  % Width
brd = 0.4;  % Breadth
Ar = wid*brd;  % Area
Iy = wid^3*brd/12;  % 2nd moment of area

%coefficients
Cb = sqrt(Ey*Iy / (rho*Ar)); %bending stiffness
Cs = sqrt(kappa*G*Ar/(rho*Ar)); %shear stiffness
Cr = sqrt(Iy / Ar); %rotational effects

% Klib K1, K2 & K3
 Klib = [struct('K', @(w,xi) sqrt( ...
     0.5*((1/Cs)^2 + (Cr/Cb)^2)*w.^2 + ...
     sqrt((w.^2)/(Cb^2) + 0.25*((1/Cs)^2 - (Cr/Cb)^2)^2 * w.^4) ));
     struct('K', @(w,xi) sqrt( ...
     (-0.5*((1/Cs)^2 + (Cr/Cb)^2)*w.^2 + ...
     sqrt((w.^2)/(Cb^2) + 0.25*((1/Cs)^2 - (Cr/Cb)^2)^2 * w.^4))));
     struct('K', @(w,xi) w*sqrt(Ey/rho))];

%Wave components
wcomps = [-1j 1; % First component -> exp(-ik1 x )
           -1 2; % Second component-> exp(-k2 x )
           1j 1; % Third component -> exp(ik1 x )
            1 2; % Fourth component-> exp(k2 x )
          -1j 3; % Fifth component-> exp(-ik3 x )
          1j 3]; % Sixth component-> exp(ik3 x )


% The relations between the coefficients of deflection wave components & bending slope wave components
P = @(w,xi) (Klib(1).K(w,xi).^2 * Cs^2 - w.^2) ./ (Klib(1).K(w,xi) * Cs^2);
N = @(w,xi) (Klib(2).K(w,xi).^2 * Cs^2 + w.^2) ./ (Klib(2).K(w,xi) * Cs^2);

L = 4;  
H = 3;   

pcs = [struct('coords',[0 0 0;0 H 0],'wcomps', wcomps);
    struct('coords',[0 H 0;L H 0],'wcomps', wcomps);   
    struct('coords',[L H 0;L 0 0],'wcomps', wcomps)];  

% Setup Boundary Conditions. Fix-Fix used here.
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1 0 0; -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0; 0 0 0 0 1 1]);
    struct('i', 6, 'cofs', @(w,xi) [1 1 1 1 0 0; -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0; 0 0 0 0 1 1])];

%% Setup the Joint
cofs = @(w,xi)[-1j*P(w,xi), -1*N(w,xi), 1j*P(w,xi), 1*N(w,xi), 0, 0, 1j*P(w,xi), 1*N(w,xi), -1j*P(w,xi), -1*N(w,xi), 0, 0;
    1-((brd/2)*(-1j*P(w,xi))), 1-((brd/2)*(-1*N(w,xi))), 1-((brd/2)*(1j*P(w,xi))), 1-((brd/2)*(1*N(w,xi))), 0, 0, 0, 0, 0, 0, -1, -1;
    0, 0, 0, 0, -1, -1, 1-((brd/2)*(-1j*P(w,xi))), 1-((brd/2)*(-1*N(w,xi))), 1-((brd/2)*(1j*P(w,xi))), 1-((brd/2)*(1*N(w,xi))), 0, 0;
    -1j*(G*Ar*kappa)*(Klib(1).K(w,xi)-P(w,xi)), 1*(G*Ar*kappa)*(-Klib(2).K(w,xi)+N(w,xi)), 1j*(G*Ar*kappa)*(Klib(1).K(w,xi)-P(w,xi)), 1*(G*Ar*kappa)*(Klib(2).K(w,xi)-N(w,xi)), 0, 0, 0, 0, 0, 0, -1j*Ey*Ar*Klib(3).K(w,xi), 1j*Ey*Ar*Klib(3).K(w,xi);
    0, 0, 0, 0, 1j*Ey*Ar*Klib(3).K(w,xi), -1j*Ey*Ar*Klib(3).K(w,xi), -1j*G*Ar*kappa*(Klib(1).K(w,xi)-P(w,xi)), 1*G*Ar*kappa*(-Klib(2).K(w,xi)+N(w,xi)), 1j*G*Ar*kappa*(Klib(1).K(w,xi)-P(w,xi)), 1*G*Ar*kappa*(Klib(2).K(w,xi)-N(w,xi)), 0, 0; 
    -Ey*Iy*[-1*(Klib(1).K(w,xi)*P(w,xi)), 1*(Klib(2).K(w,xi)*N(w,xi)), -1*(Klib(1).K(w,xi)*P(w,xi)), 1*(Klib(2).K(w,xi)*N(w,xi)), 0, 0, -1*(Klib(1).K(w,xi)*P(w,xi)), 1*(Klib(2).K(w,xi)*N(w,xi)), -1*(Klib(1).K(w,xi)*P(w,xi)), 1*(Klib(2).K(w,xi)*N(w,xi)), 0, 0] +...
    (brd*0.5*G*Ar*kappa)*[-1j*(Klib(1).K(w,xi)-P(w,xi)), 1*(-Klib(2).K(w,xi)+N(w,xi)), 1j*(Klib(1).K(w,xi)-P(w,xi)), 1*(Klib(2).K(w,xi)-N(w,xi)), 0, 0, -1j*(Klib(1).K(w,xi)-P(w,xi)), 1*(-Klib(2).K(w,xi)+N(w,xi)), 1j*(Klib(1).K(w,xi)-P(w,xi)), 1*(Klib(2).K(w,xi)-N(w,xi)), 0, 0]
    ];

joints = [struct('type', 2, 'i', 2, 'j', 3, 'cofs', cofs);
          struct('type', 2, 'i', 4, 'j', 5, 'cofs',cofs)];
    
%%
[pcs, bcs, joints, ~, Klib] = WBPREPROC(pcs, bcs, joints, [], Klib);

%% 
Nw = 500;
Ws = linspace(0, 1e1, Nw+1);
Ws = Ws(2:end);
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