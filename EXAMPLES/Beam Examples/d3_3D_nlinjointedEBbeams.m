clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a nonlinear jointed 3D/spatial Euler-Bernoulli Beams example
%Parameters are taken from Ljiljana, et al, 2016.

%% Setup Model
Ey = 3e10;
nu = 0.25;
G = Ey/(2*(1+nu));
rho = 2548.42;
wid = 0.3; 
brd = 0.4;  
Ar = wid*brd;  
Iy = wid^3*brd/12; 
Iz = brd^3*wid/12;  
L0 = 8.0; 
Klib = [struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25)); %wavenumber for bending in y direction
        struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iz)^(0.25)); %wavenumber for bending in z direction
        struct('K', @(w,xi) w*sqrt(rho/Ey))]; %torisonal wavenumber

wcomps = [1 1;  % exp(k1 x )
         -1 1;  % exp(-k1 x )
         1j 1;  % exp(ik x )
        -1j 1;  % exp(ik1 x)
          1 2;  % exp(k2 x )
         -1 2;  % exp(-k2 x )
         1j 2;  % exp(ik2 x )
        -1j 2;  % exp(-ik2 x )
         1j 3;  % exp(ik3 x )
        -1j 3]; % exp(-ik3 x )

% Setup "wave-based pieces"
pcs = [struct('coords', [0 0 0;L0/3 0 0;L0/2 0 0], 'wcomps', wcomps);
    struct('coords', [L0/2 0 0; L0 0 0], 'wcomps', wcomps)];
% Setup Boundary Conditions. Fix-Fix used here.
bcs = [struct('i', 1, 'cofs', @(w,xi) [0 0 0 0 0 0 0 0 1 1; 1 1 1 1 0 0 0 0 0 0; 0 0 0 0 1 1 1 1 0 0;1 -1 1j -1j 0 0 0 0 0 0 ;0 0 0 0 1 -1 1j -1j 0 0]);
    struct('i', 5, 'cofs', @(w,xi) [0 0 0 0 0 0 0 0 1 1; 1 1 1 1 0 0 0 0 0 0; 0 0 0 0 1 1 1 1 0 0;1 -1 1j -1j 0 0 0 0 0 0 ;0 0 0 0 1 -1 1j -1j 0 0])];

%% Excitation
Mx = @(w,xi)inv([(-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0 0 0 0 0];
                 (-Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0 0 0 0 0];
                  [1 1 1 1 0 0 0 0 0 0];
                  Klib(1).K(w,xi)*[1 -1 1j -1j 0 0 0 0 0 0];
                  (-Ey*Iz*Klib(2).K(w,xi)^3)*[0 0 0 0 1 -1 -1j 1j 0 0];
                  (-Ey*Iz*Klib(2).K(w,xi)^2)*[0 0 0 0 1 1 -1 -1 0 0];
                  [0 0 0 0 1 1 1 1 0 0];
                  Klib(2).K(w,xi)*[0 0 0 0 1 -1 1j -1j 0 0];
                  (Ey*Ar*Klib(3).K(w,xi))*[0 0 0 0 0 0 0 0 1j -1j];
                  [0 0 0 0 0 0 0 0 1 1]]);
excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0;1/2;0;0;0;0;0], ...
    'rcofs0', [1/2;0;0;0;1/2;0;0;0;0;0]);
%% Setup Joint
h = [1; 3];
Nt = 128;
kJs = diag([1e9 1e9 1e9 1e9 1e9]);
cJs = diag([9e4 9e4 9e4 9e4 9e4]);
gJs = diag([2e12 2e12 2e12 2e12 0]);
cofs = @(w,xi) [(-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    (-Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    (-Ey*Iz*Klib(2).K(w,xi)^3)*[0 0 0 0 1 -1 -1j 1j 0 0 0 0 0 0 0 0 0 0 0 0];
    (-Ey*Iz*Klib(2).K(w,xi)^2)*[0 0 0 0 1 1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0];
    (Ey*Ar*Klib(3).K(w,xi))*[0 0 0 0 0 0 0 0 1j -1j 0 0 0 0 0 0 0 0 0 0];
    (-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0 0 0 0 0 -1 1 1j -1j 0 0 0 0 0 0];
    (-Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0 0 0 0 0 -1 -1 1 1 0 0 0 0 0 0];
    (-Ey*Iz*Klib(2).K(w,xi)^3)*[0 0 0 0 1 -1 -1j 1j 0 0 0 0 0 0 -1 1 1j -1j 0 0];
    (-Ey*Iz*Klib(2).K(w,xi)^2)*[0 0 0 0 1 1 -1 -1 0 0 0 0 0 0 -1 -1 1 1 0 0];
    (Ey*Ar*Klib(3).K(w,xi))*[0 0 0 0 0 0 0 0 1j -1j 0 0 0 0 0 0 0 0 -1j 1j]];
joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJs, cJs, gJs, h, Nt), ...
    'nldcofs', @(w,xi) [1 1 1 1 0 0 0 0 0 0 -1 -1 -1 -1 0 0 0 0 0 0; Klib(1).K(w,xi)*[1 -1 1j -1j 0 0 0 0 0 0 -1 1 -1j 1j 0 0 0 0 0 0]; ...
    0 0 0 0 1 1 1 1 0 0 0 0 0 0 -1 -1 -1 -1 0 0; Klib(2).K(w,xi)*[0 0 0 0 1 -1 1j -1j 0 0 0 0 0 0 -1 1 -1j 1j 0 0];0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 -1 1], ...
    'nlfcofs', @(w,xi) [eye(5);zeros(5)]);
%NOTE: This is a joint that engages both the transverse displacement as
%well as rotation such that,
%       fnl = kJs [u;th] + cJs [udot;thdot] + gJs [u^3;th^3];
%   represents the nonlinear force. "th" in the above is the rotation. 

%% Pre-Processing
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%% Conduct Nonlinear Forced Response Analysis using Continuation
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);

Wst = 138.6;
Wen = 139.4;
dw = 0.8;

Copt = struct('Nmax', 700, 'angopt', 1e-2, 'DynDscale', 1);
Famps = 2e2*[1 5 10];
acC = cell(size(Famps));
for fi=1:length(Famps)
    % Setup Linear Initial Guess
    %[Amat, ~, ~, Fv] = WVAMAT([Wst;0], h, pcs, bcs, joints, Klib, 'r');
    %ari0 = Amat\Fv*Famps(fi); 
    %NOTE: This does NOT linearize the joint. This merely assumes no joint
    %is present.
    ari0 = zeros(Npts*Nwc*Nhc, 1);
    Copt.Dscale = [abs(ari0+1e-5);Wst];  % Setup Scaling for continuation (sometimes helps)

    % Conduct continuation
    ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famps(fi), h, pcs, bcs, joints, Klib), ...
        ari0, Wst, Wen, dw, Copt);

    % Convert to complex representation
    acC{fi} = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
    acC{fi}([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];
end

%% Plot Results
opi = 11:20;
figure(1)
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

