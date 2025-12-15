clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a nonlinear jointed Timoshenko Beam example

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
L0 = 2.0;  % Total Length

%coefficients
Cb = sqrt(Ey*Iy / (rho*Ar)); %bending stiffness
Cs = sqrt(kappa*G*Ar/(rho*Ar)); %shear stiffness
Cr = sqrt(Iy / Ar); %rotational effects

%Wavenumbers K1 & K2
 Klib = [struct('K', @(w,xi) sqrt( ...
     0.5*((1/Cs)^2 + (Cr/Cb)^2)*w.^2 + ...
     sqrt((w.^2)/(Cb^2) + 0.25*((1/Cs)^2 - (Cr/Cb)^2)^2 * w.^4) ));
     struct('K', @(w,xi) sqrt( ...
     abs(0.5*((1/Cs)^2 + (Cr/Cb)^2)*w.^2 - ...
     sqrt((w.^2)/(Cb^2) + 0.25*((1/Cs)^2 - (Cr/Cb)^2)^2 * w.^4))))];

%Wave components
wcomps = [-1j 1;  % First component -> exp(-ik1 x )
         -1 2;  % Second component-> exp(-k2 x )
         1j 1;  % Third component -> exp(ik1 x )
        1 2]; % Fourth component-> exp(k2 x )

%pieces
pcs = [struct('coords', [0;L0/3;L0], 'wcomps', wcomps);
    struct('coords', [L0;2*L0], 'wcomps', wcomps)];

% The relations between the coefficients of deflection wave components & bending slope wave components
P = @(w, xi) Klib(1).K(w, xi) * (1 - ((w.^2) / (Klib(1).K(w, xi).^2 * Cs^2)));
N = @(w, xi) Klib(2).K(w, xi) * (1 + ((w.^2) / (Klib(2).K(w, xi).^2 * Cs^2)));

% Fixed-Fixed Boundary conditions 
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1; -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi)]);
    struct('i', 5, 'cofs', @(w,xi) [1 1 1 1; -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi)])];

%% Excitation
Mx = @(w,xi)inv([G*Ar*kappa*[1j*(-Klib(1).K(w,xi)+P(w,xi)) 1*(-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) 1*(Klib(2).K(w,xi)-N(w,xi))];
                  -Ey*Iy*[-1*(Klib(1).K(w,xi).*P(w,xi)) 1*(Klib(2).K(w,xi).*N(w,xi)) -1*(Klib(1).K(w,xi).*P(w,xi)) 1*(Klib(2).K(w,xi)*N(w,xi))];
                  [1 1 1 1];
                  [-1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi)]]);

excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0], ...
    'rcofs0', [1/2;0;0;0]);

%% joint setup
h = [1; 3];
Nt = 128;

kJs = diag([1e9 1e9]);
cJs = diag([320 320]);
gJs = diag([1e8 0]);

cofs = @(w,xi)[G*Ar*kappa*[1j*(-Klib(1).K(w,xi)+P(w,xi)) 1*(-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) 1*(Klib(2).K(w,xi)-N(w,xi)) 0 0 0 0];
    -Ey*Iy*[-1*(Klib(1).K(w,xi)*P(w,xi)) 1*(Klib(2).K(w,xi)*N(w,xi)) -1*(Klib(1).K(w,xi)*P(w,xi)) 1*(Klib(2).K(w,xi)*N(w,xi)) 0 0 0 0];
    G*Ar*kappa*[1j*(-Klib(1).K(w,xi)+P(w,xi)) 1*(-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) 1*(Klib(2).K(w,xi)-N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) 1*(Klib(2).K(w,xi)-N(w,xi)) 1j*(-Klib(1).K(w,xi)+P(w,xi)) 1*(-Klib(2).K(w,xi)+N(w,xi))];
    -Ey*Iy*[-1*(Klib(1).K(w,xi)*P(w,xi)) 1*(Klib(2).K(w,xi)*N(w,xi)) -1*(Klib(1).K(w,xi)*P(w,xi)) 1*(Klib(2).K(w,xi)*N(w,xi)) 1*(Klib(1).K(w,xi)*P(w,xi)) -1*(Klib(2).K(w,xi)*N(w,xi)) 1*(Klib(1).K(w,xi)*P(w,xi)) -1*(Klib(2).K(w,xi)*N(w,xi))]];

joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJs, cJs, gJs, h, Nt), ...
    'nldcofs', @(w,xi) [1 1 1 1 -1 -1 -1 -1; -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) -1j*P(w,xi) -1*N(w,xi)], ...
    'nlfcofs', @(w,xi) [eye(2);zeros(2)]);

%% Pre-Processing
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%% Conduct Nonlinear Forced Response Analysis using Continuation
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);

Wst = 1020.8; 
Wen = 1021.2;
dw = 0.5;

Copt = struct('Nmax', 300, 'angopt', 1e-1, 'DynDscale', 1);
Famps = 2e3*[1 10 20];
acC = cell(size(Famps));
for fi=1:length(Famps)
    % Setup Linear Initial Guess
    [Amat, ~, ~, Fv] = WVAMAT([Wst;0], h, pcs, bcs, joints, Klib, 'r');
    ari0 = Amat\Fv*Famps(fi); 
    %ari0 = zeros(Npts*Nwc*Nhc, 1);
    Copt.Dscale = [abs(ari0+1e-6);Wst];  % Setup Scaling for continuation (sometimes helps)
    % Conduct continuation
    ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famps(fi), h, pcs, bcs, joints, Klib), ...
        ari0, Wst, Wen, dw, Copt);
    % Convert to complex representation
    acC{fi} = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
    acC{fi}([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];
end

%% Plot Results
opi = 5:8;
figure(2)
clf()
aa = gobjects(size(Famps));
for fi=1:length(Famps)
    subplot(2,1,1)
    aa(fi)=plot(acC{fi}(end,:), abs(sum(2*acC{fi}(opi,:))), '-', 'LineWidth', 2); hold on
    legend(aa(fi), sprintf('F = %.0f kN', Famps(fi)/1e3));
    grid on
    ylabel('Response (m)')
    subplot(2,1,2)
    plot(acC{fi}(end,:), rad2deg(angle(sum(2*acC{fi}(opi,:)))), '-', 'LineWidth', 2); hold on
    grid on
    ylabel('Phase (degs)')
end
subplot(2,1,1)
xlim([Wst Wen])
legend(aa, 'Location', 'northwest')
subplot(2,1,2)
xlim([Wst Wen])
set(gca, 'YTick', -180:90:180)
xlabel('Frequency (rad/s)')
