clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')
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
%pieces
L = 4.0;  
H = 3.0;
pcs = [struct('coords',[0 0;0 2.5;0 H],'wcomps', wcomps);
    struct('coords',[0 H;L H],'wcomps', wcomps);   
    struct('coords',[L H;L 0],'wcomps', wcomps)]; 
% BCs (Fix_Fix)
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1 0 0; -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0; 0 0 0 0 1 1]);
    struct('i', 7, 'cofs', @(w,xi) [1 1 1 1 0 0; -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0; 0 0 0 0 1 1])];
%joint
cofs = @(w,xi)[1+((brd/2)*(-1j*P(w,xi))), 1+((brd/2)*(-1*N(w,xi))), 1+((brd/2)*(1j*P(w,xi))), 1+((brd/2)*(1*N(w,xi))), 0, 0, 0, 0, 0, 0, 1, 1;
    0, 0, 0, 0, -1, -1, 1+((brd/2)*(-1j*P(w,xi))), 1+((brd/2)*(-1*N(w,xi))), 1+((brd/2)*(1j*P(w,xi))), 1+((brd/2)*(1*N(w,xi))), 0, 0;
    -1j*P(w,xi), -1*N(w,xi), 1j*P(w,xi), 1*N(w,xi), 0, 0, 1j*P(w,xi), 1*N(w,xi), -1j*P(w,xi), -1*N(w,xi), 0, 0;
    -1j*(G*Ar*kappa)*(Klib(1).K(w,xi)-P(w,xi)), 1*(G*Ar*kappa)*(-Klib(2).K(w,xi)+N(w,xi)), 1j*(G*Ar*kappa)*(Klib(1).K(w,xi)-P(w,xi)), 1*(G*Ar*kappa)*(Klib(2).K(w,xi)-N(w,xi)), 0, 0, 0, 0, 0, 0, -1j*Ey*Ar*Klib(3).K(w,xi), 1j*Ey*Ar*Klib(3).K(w,xi);
    0, 0, 0, 0, 1j*Ey*Ar*Klib(3).K(w,xi), -1j*Ey*Ar*Klib(3).K(w,xi), -1j*(G*Ar*kappa)*(Klib(1).K(w,xi)-P(w,xi)), 1*(G*Ar*kappa)*(-Klib(2).K(w,xi)+N(w,xi)), 1j*(G*Ar*kappa)*(Klib(1).K(w,xi)-P(w,xi)), 1*(G*Ar*kappa)*(Klib(2).K(w,xi)-N(w,xi)), 0, 0; 
    -Ey*Iy*[-1*(Klib(1).K(w,xi)*P(w,xi)), 1*(Klib(2).K(w,xi)*N(w,xi)), -1*(Klib(1).K(w,xi)*P(w,xi)), 1*(Klib(2).K(w,xi)*N(w,xi)), 0, 0, 1*(Klib(1).K(w,xi)*P(w,xi)), -1*(Klib(2).K(w,xi)*N(w,xi)), 1*(Klib(1).K(w,xi)*P(w,xi)), -1*(Klib(2).K(w,xi)*N(w,xi)), 0, 0] +...
    (brd*0.5*G*Ar*kappa)*[-1j*(Klib(1).K(w,xi)-P(w,xi)), 1*(-Klib(2).K(w,xi)+N(w,xi)), 1j*(Klib(1).K(w,xi)-P(w,xi)), 1*(Klib(2).K(w,xi)-N(w,xi)), 0, 0, -1j*(Klib(1).K(w,xi)-P(w,xi)), 1*(-Klib(2).K(w,xi)+N(w,xi)), 1j*(Klib(1).K(w,xi)-P(w,xi)), 1*(Klib(2).K(w,xi)-N(w,xi)), 0, 0]
    ];

joints = [struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs);
          struct('type', 2, 'i', 5, 'j', 6, 'cofs', cofs)];

%% Excitation
Mx = @(w,xi)inv([G*Ar*kappa*[1j*(-Klib(1).K(w,xi)+P(w,xi)) 1*(-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) 1*(Klib(2).K(w,xi)-N(w,xi)) 0 0];
                  -Ey*Iy*[-1*(Klib(1).K(w,xi).*P(w,xi)) 1*(Klib(2).K(w,xi).*N(w,xi)) -1*(Klib(1).K(w,xi).*P(w,xi)) 1*(Klib(2).K(w,xi)*N(w,xi)) 0 0];
                  Ey*Ar*Klib(3).K(w,xi)*[0 0 0 0 -1j 1j];
                  [1 1 1 1 0 0];
                  [-1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0];
                  [0 0 0 0 1 1]]);

excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0;0;0], 'rcofs0', [1/2;0;0;0;0;0]);

%% Preprocess Everything
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);
Nwc = size(wcomps,1);  % Number of wave components
%% Conduct Linear Forced Response Analysis
Nw = 500;
Ws = linspace(1, 1e3, Nw);
Npts = pcs(end).irange(end);  
ACs = zeros(Npts*Nwc,Nw);
Famp = 200;
for iw=1:Nw
    [Amat, ~, ~, Fv] = WVAMAT([Ws(iw);0], 1, pcs, bcs, joints, Klib);
    ACs(:,iw) = Amat\(Fv*Famp);
end
%% Plot Forced Response
opi = 7:12;  % Output wave coefficients
figure(1)
clf()
subplot(2,1,1)
semilogy(Ws/1e3, abs(2*sum(ACs(opi,:))));
grid on
ylabel('Response (m)')
subplot(2,1,2)
plot(Ws/1e3, rad2deg((angle(2*sum(ACs(opi,:))))));
grid on
set(gca, 'YTick', -180:90:180)
ylabel('Phase (degs)')
xlabel('Frequency (k rad/s)')