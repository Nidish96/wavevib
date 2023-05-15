clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a nonlinear jointed Euler-Bernoulli Beam example
%similar to the jointed bar example. Parameters taken from 
% Krishna and Chandramouli, 2012

%% Setup Model
Ey = 2.1e11;
rho = 7680;
thk = 3e-3;   % Thickness 1
wid  = 40e-3;  % Width
Ar = thk*wid;  % Area
Iy = thk^3*wid/12;  % 2nd moment of area
L1  = 560e-3;  % Primary beam length
xF1 = L1/6;    % Excitation Location (on primary beam)
L2 = 445e-3;  % Secondary beam length (each)
% Klib = struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25));  % Undamped EB-Beam Case
al = 0.80;    % 0.80, 1.80
bt = 1.1e-4;  % 1.1e-4, 2.475e-4
Klib = struct('K', @(w,xi) ((rho*Ar*(w.^2+1j*w*al))./(Ey*Iy*(1-1j*w*bt))).^(0.25) );
wcomps = [1 1;  % First component -> exp(  k x )
         -1 1;  % Second component-> exp( -k x )
         1j 1;  % Third component -> exp( ik x )
        -1j 1]; % Fourth component-> exp(-ik x )

% Setup "wave-based pieces"
pcs = [struct('coords', [0;xF1;L1], 'wcomps', wcomps);
    struct('coords', [L1;L1+L2], 'wcomps', wcomps);
    struct('coords', [L1;L1+L2], 'wcomps', wcomps)];

% Setup Boundary Conditions. Fix-Fix used here.
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j]);  % fixed end
    struct('i', 5, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j]);     % fixed end
    struct('i', 7, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j]);     % fixed end
    struct('i', 3, 'cofs', @(w,xi) [1 1 -1 -1]);        % Moment-free
    struct('i', 4, 'cofs', @(w,xi) [1 1 -1 -1]);        % Moment-free
    struct('i', 6, 'cofs', @(w,xi) [1 1 -1 -1])];       % Moment-free

%% Setup Excitation
Mx = @(w,xi) inv([Ey*Iy*Klib.K(w,xi)^3*[-1 1 1j -1j];
                  Ey*Iy*Klib.K(w,xi)^2*[-1 -1 1 1];
                  [1 1 1 1];
                  Klib.K(w,xi)*[1 -1 1j -1j]]);

excs = struct('i', 2, 'nh', 1, ...
    'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0]);
%'nh' sets the harmonic at which to apply the excitation

%% Setup AFT parameters
h = (0:5)';
Nt = 2^10;

%% Setup Joints
knl = 220/2;      % 220,      880,    (242, 484, 880, 1210)
gap = 2.5e-3;   % 2.5e-3,   0.35,   0.35

cofs = @(w,xi) [-[zeros(1,4) 1 -1 -1j 1j zeros(1,4)];
    -[zeros(1,8) 1 -1 -1j 1j];
    kron([1 -1 -1], [1 -1 -1j 1j])];
joints = struct('type', 3, 'is', [3 4 6], ...
    'cofs', cofs, 'nl', @(Uw) HCONTACT(Uw, knl, gap, h, Nt), ...
    'nldcofs', @(w,xi) kron([1 -1 0;1 0 -1], [1 1 1 1]), ...
    'nlfcofs', @(w,xi) [eye(2);0 0]/Klib.K(w,xi)^3);

%% Pre-Processing
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%% Conduct Nonlinear Forced Response Analysis using Continuation
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);

Wst = 2*pi*4;
Wen = 2*pi*12;
dw = 0.25;

opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');

Copt = struct('Nmax', 600, 'angopt', 1e-1, 'DynDscale', 1);
Famps = [0.4 1.0 2.5];  % 0.336, 0.336
acC = cell(size(Famps));
for fi=1:length(Famps)
    [Amat, ~, ~, Fv, ~, ~, JEV] = WVAMAT([Wst;0], h, pcs, bcs, joints, Klib, 'r');
    return;
    ari0 = Amat\Fv;

    Copt.Dscale = [1e-6*ones(size(ari0));Wst];
    ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famps(fi), h, pcs, bcs, joints, Klib), ...
        ari0, Wst, Wen, dw, Copt);

    % Convert to complex representation
    acC{fi} = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
    acC{fi}([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];
end
%% Plot Results
opi = Npts*4*(h(1)==0)+(13:16);
figure(1)
clf()
aa = gobjects(size(Famps));
for fi=1:length(Famps)
    subplot(2,1,1)
    aa(fi)=plot(acC{fi}(end,:)/2/pi, abs(sum(2*acC{fi}(opi,:))), '.-'); hold on
    legend(aa(fi), sprintf('F = %f N', Famps(fi)));
    grid on
    ylabel('Response (m)')
    subplot(2,1,2)
    plot(acC{fi}(end,:)/2/pi, rad2deg(angle(sum(2*acC{fi}(opi,:)))), '.-'); hold on
    grid on
    ylabel('Phase (degs)')
end
subplot(2,1,1)
xlim(sort([Wst Wen]/2/pi))
legend(aa, 'Location', 'northwest')
subplot(2,1,2)
xlim(sort([Wst Wen]/2/pi))
set(gca, 'YTick', -180:90:180)
ylim([-1 1]*180)
xlabel('Frequency (rad/s)')
