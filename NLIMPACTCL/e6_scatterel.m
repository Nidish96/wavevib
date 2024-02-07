clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')
addpath('../ROUTINES/FOCHLERD/')

addpath('../ROUTINES/export_fig')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a nonlinear jointed Euler-Bernoulli Beam example
%similar to the jointed bar example. Parameters taken from 
% Krishna and Chandramouli, 2012

savfig = false;
animfig = false;
savdat = false;
analyze = true;
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
              % Joint parameters
knl = 1210/2;      % 220,      880,    (242, 484, 880, 1210)
gap = 2.5e-3;   % 2.5e-3,   0.35,   0.35

Klib = struct('K', @(w,xi) ((rho*Ar*(w.^2+1j*w*al))./(Ey*Iy*(1-1j*w*bt))).^(0.25) );
wcomps = [1 1;  % First component -> exp(  k x )
          -1 1;  % Second component-> exp( -k x )
          1j 1;  % Third component -> exp( ik x )
          -1j 1]; % Fourth component-> exp(-ik x )

% Setup "wave-based pieces"
pcs = [struct('coords', [0 0;xF1 0;L1 0], 'wcomps', wcomps);
       struct('coords', [L1 gap+thk;L1+L2 gap+thk], 'wcomps', wcomps);
       struct('coords', [L1 -gap-thk;L1+L2 -gap-thk], 'wcomps', wcomps)];

% Setup Boundary Conditions. Fix-Fix used here.
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j], ...
              'cofs0', @(xi) [eye(2) zeros(2)]);  % fixed end
       struct('i', 5, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j], ...
              'cofs0', @(xi) [eye(2) zeros(2)]);     % fixed end
       struct('i', 7, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j], ...
              'cofs0', @(xi) [eye(2) zeros(2)]);     % fixed end
       struct('i', 3, 'cofs', @(w,xi) [1 1 -1 -1], ...
              'cofs0', @(xi) [0 0 1 0]);        % Moment-free
       struct('i', 4, 'cofs', @(w,xi) [1 1 -1 -1], ...
              'cofs0', @(xi) [0 0 1 0]);        % Moment-free
       struct('i', 6, 'cofs', @(w,xi) [1 1 -1 -1], ...
              'cofs0', @(xi) [0 0 1 0])];       % Moment-free


%% Setup Excitation
Mx = @(w,xi) inv([Ey*Iy*Klib.K(w,xi)^3*[-1 1 1j -1j];
                  Ey*Iy*Klib.K(w,xi)^2*[-1 -1 1 1];
                  [1 1 1 1];
                  Klib.K(w,xi)*[1 -1 1j -1j]]);
excs = struct('i', 2, 'nh', 1, ...
              'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0], ...
              'rcofs0', @(xi) [0;0;0;1/(6*Ey*Iy)]);
%'nh' sets the harmonic at which to apply the excitation

%% Setup AFT parameters
h = (0:16)'; % (1:16);
Nt = 2^10;

%% Setup Joints
cofs = @(w,xi) [-[zeros(1,4) 1 -1 -1j 1j zeros(1,4)];
                -[zeros(1,8) 1 -1 -1j 1j];
                kron([1 -1 -1], [1 -1 -1j 1j])];
cofs0 = @(xi) [-[zeros(1,4) 0 0 0 6 zeros(1,4)];
               -[zeros(1,8) 0 0 0 6];
               kron([1 -1 -1], [0 0 0 1])];
joints = struct('type', 3, 'is', [3 4 6], ...
                'cofs', cofs, 'nl', @(Uw) HCONTACT(Uw, knl, gap, h, Nt), ...
                'nldcofs', @(w,xi) kron([1 -1 0;-1 0 1], [1 1 1 1]), ...
                'nlfcofs', @(w,xi) [1 0;0 -1;0 0]/(Ey*Iy*Klib.K(w,xi)^3), ...
                ...
                'cofs0', cofs0, ...
                'nldcofs0', @(xi) kron([1 -1 0;-1 0 1], [1 0 0 0]), ...
                'nlfcofs0', @(xi) [1 0;0 -1;0 0]/(Ey*Iy));

%% Pre-Processing
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);
% Klib.dKdw = @(w,xi) Klib.dKdw(w+eps,xi);  % Adjusting for w=0.

%% Load Forced Response Results
Npcs = length(pcs);
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);

load(sprintf('./DATS/E_NLIMPDAT_2Knl%d_H%d.mat', knl,Nh), ...
     'acC', 'pds', 'h', 'Famps', 'Klib', 'wcomps', ...
     'bcs', 'joints', 'excs', ...
     'knl', 'gap', 'Wst', 'Wen')
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);
Nf = length(Famps);


stabis = cellfun(@(pp) sign(conv(sign(pp(:)')==sign(pp(1)), [1 1], 'same')), ...
                 pds, 'Uniformoutput', false);
opis = {(13:16), (17:20), (25:28)};

%% Estimations of joint scattering

