clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a EB beam with a nonlinear single lap bolted joint 
% under various forcing example

%% Setup
Ey = 210e9;
rho = 7850;
wid = 1e-2;  % Width
brd = 1e-2;  % Breadth
Ar = wid*brd;  % Area
Iy = wid^3*brd/12;  % 2nd moment of area
L0 = 0.5;  % Total Length
alpha = 0.25; %mass proportional damping coeff
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

%% Nonlinear Joint setup
h = [1;3]; % List of harmonics 
Nt = 128; % Number of AFT points
mu = 0.2; %friction coefficient
gap = -1e-6; %interface's initial gap

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

%tangential, normal, and rotational stiffness 
kt = 1/((((h1+h2)/(2*dia))^a) * (b/n) * ((1/(h1*E1)) + (1/(n*h2*E2)) ...
    + (1/(n*h1*Eb)) + (1/(2*n*h2*Eb)))); %Huth's Formula
k_bolt = (Eb * pi * (dia/2)^2) /(h1+h2);
k_comp = 0.5774*pi*Eb*dia/(2*log(5*((0.5774*(h1+h2)) ...
    +(0.5*dia))/((0.5774*(h1+h2))+(2.5*dia))));
kn = k_comp + k_bolt;
k_theta = kn*lj^2 / 12;

%Spring-slider Interaction forces and moments, Joint Equilibiriums
cofs = @(w,xi)[Ey*Ar*Klib(2).K(w,xi)*[0 0 0 0 1j -1j 0 0 0 0 0 0];  % Axial force left = Ft
    (-Ey*Iy*Klib(1).K(w,xi)^3)*[-1j 1j 1 -1 0 0 0 0 0 0 0 0];  % Shear force left = Fn
    (-Ey*Iy*Klib(1).K(w,xi)^2)*[0 0 0 0 0 0 -1 -1 1 1 0 0] ...
    - k_theta*[-1j 1j -1 1 0 0 1j -1j 1 -1 0 0]; % Moment left = Fn x e
    Ey*Ar*Klib(2).K(w,xi)*[0 0 0 0 1j -1j 0 0 0 0 -1j 1j];  % Axial force left-right=0
    (-Ey*Iy*Klib(1).K(w,xi)^3)*[-1j 1j 1 -1 0 0 1j -1j -1 1 0 0];  % Shear force left-right=0
    (-Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0 -1 -1 1 1 0 0] %M2-M1-h/2(F1+F2)+d/2(V1+V2)=0
    ];

joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) HFRIC(Uw, mu, kt, kn, gap, h, Nt), ...
    'nldcofs', @(w,xi) [[0 0 0 0 1 1 0 0 0 0 -1 -1]-(0.5*brd*[1j -1j 1 -1 0 0 1j -1j 1 -1 0 0]); ...
    [1 1 1 1 0 0 -1 -1 -1 -1 0 0]], ...
    'nlfcofs', @(w,xi) [1 0;0 1;(-0.5*h1) 0;0 0;0 0;h1 -dia]);
%nldcofs: (u1-u2) - h/2(psi1+psi2); y1-y2
%nlfcofs: Ft 0; 0 Fn; (-h/2)*Ft 0;0 0; 0 0;h*Ft -d*Fn

%% Excitation
Mx = @(w,xi)inv([(-Ey*Iy*Klib(1).K(w,xi)^3)*[-1j 1j 1 -1 0 0];
                 (-Ey*Iy*Klib(1).K(w,xi)^2)*[-1 -1 1 1 0 0];
                 (Ey*Ar*Klib(2).K(w,xi))*[0 0 0 0 1j -1j];
                  [1 1 1 1 0 0];
                  [1j -1j 1 -1 0 0];
                  [0 0 0 0 1 1]]);

excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0;0;0], ...
    'rcofs0', [1/2;0;0;0;0;0]);

%% Pre-Processing
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%% Conduct Nonlinear Forced Response Analysis using Continuation
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);

%Freq. Sweeping 
Wen = 3400;  
Wst = 4000; 
dw = 0.05;

Copt = struct('Nmax', 5000, 'angopt', 1e-1, 'DynDscale', 1, 'solverchoice', 3);
Famps = [1 10 25 50 100];
acC = cell(size(Famps));
for fi=1:length(Famps)
    % Setup Linear Initial Guess
    %ari0 = zeros(Npts*Nwc*Nhc, 1);
    [Amat, ~, ~, Fv] = WVAMAT([Wst;0], h, pcs, bcs, joints, Klib, 'r');
    ari0 = Amat\Fv*Famps(fi); 
    %NOTE: This does NOT linearize the joint. This merely assumes no joint
    %is present.
    Copt.Dscale = [abs(ari0+1e-5);Wst];  % Setup Scaling for continuation (sometimes helps)

    % Conduct continuation
    ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famps(fi), h, pcs, bcs, joints, Klib), ...
        ari0, Wst, Wen, dw, Copt);

    % Convert to complex representation
    acC{fi} = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
    acC{fi}([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];
end

%% Plot Results
opi = 7:10;
figure()
clf;
hold on
aa = gobjects(size(Famps));
for fi=1:length(Famps)
    subplot(2,1,1)
    aa(fi)=semilogy(acC{fi}(end,:), (abs(sum(2*acC{fi}(opi,:)))), '-', 'LineWidth', 2); hold on
    legend(aa(fi), sprintf('F = %.0f N', Famps(fi)));
    grid on
    ylabel('Response (m)')
    subplot(2,1,2)
    plot(acC{fi}(end,:), rad2deg(angle(sum(2*acC{fi}(opi,:)))), '-', 'LineWidth', 2); hold on
    grid on
    ylabel('Phase (degs)')
end
subplot(2,1,1)
xlim(sort([Wst Wen]))
legend(aa, 'Location', 'northwest')
subplot(2,1,2)
xlim(sort([Wst Wen]))
set(gca, 'YTick', -180:90:180)
xlabel('Frequency (rad/s)')
