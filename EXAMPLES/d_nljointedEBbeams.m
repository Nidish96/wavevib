clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a nonlinear jointed Euler-Bernoulli Beam example
%similar to the jointed bar example.

%% Setup Model
Ey = 190e9;
rho = 7680;
wid = 0.2;  % Width
brd = 0.4;  % Breadth
Ar = wid*brd;  % Area
Iy = wid^3*brd/12;  % 2nd moment of area
L0 = 2.0;  % Total Length
Klib = struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25));
wcomps = [1 1;  % First component -> exp(  k x )
         -1 1;  % Second component-> exp( -k x )
         1j 1;  % Third component -> exp( ik x )
        -1j 1]; % Fourth component-> exp(-ik x )

% Setup "wave-based pieces"
pcs = [struct('coords', [0;L0/3;L0], 'wcomps', wcomps);
    struct('coords', [L0;2*L0], 'wcomps', wcomps)];

% Setup Boundary Conditions. Fix-Fix used here.
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1]);
    struct('i', 1, 'cofs', @(w,xi) [1 -1 1j -1j]);
    struct('i', 5, 'cofs', @(w,xi) [1 1 1 1]);
    struct('i', 5, 'cofs', @(w,xi) [1 -1 1j -1j])];

%% Setup Excitation
Mx = @(w,xi) inv([Ey*Iy*Klib.K(w,xi)^3*[-1 1 1j -1j];
                  Ey*Iy*Klib.K(w,xi)^2*[-1 -1 1 1];
                  [1 1 1 1];
                  Klib.K(w,xi)*[1 -1 1j -1j]]);

excs = struct('i', 2, 'nh', 1, ...
    'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0]);
%'nh' sets the harmonic at which to apply the excitation

%% Setup Joint
h = [1; 3];
Nt = 128;

kJs = diag([1e9 1e9]);
cJs = diag([320 320]);
gJs = diag([1e8 0]);
cofs = @(w,xi) [-Klib.K(w,xi)^3*[1 -1 -1j 1j 0 0 0 0];
    -Klib.K(w,xi)^2*[1 1 -1 -1 0 0 0 0];
    -Klib.K(w,xi)^3*[1 -1 -1j 1j -1 1 1j -1j];
    -Klib.K(w,xi)^2*[1 1 -1 -1 -1 -1 1 1]];
joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJs, cJs, gJs, h, Nt), ...
    'nldcofs', @(w,xi) [1 1 1 1 -1 -1 -1 -1; Klib.K(w,xi)*[1 -1 1j -1j -1 1 -1j 1j]], ...
    'nlfcofs', @(w,xi) [eye(2);zeros(2)]/(Ey*Iy));
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

Wst = 1055.5;
Wen = 1055.9;
dw = 0.1;

Copt = struct('Nmax', 300, 'angopt', 1e-1, 'DynDscale', 1);
Famps = 2e3*[1 10 20];
acC = cell(size(Famps));
for fi=1:length(Famps)
    % Setup Linear Initial Guess
    [Amat, ~, ~, Fv] = WVAMAT([Wst;0], h, pcs, bcs, joints, Klib, 'r');
    ari0 = Amat\Fv*Famps(fi); 
    %NOTE: This does NOT linearize the joint. This merely assumes no joint
    %is present.
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
