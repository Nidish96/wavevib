clc
clear all
addpath('../ROUTINES/SOLVERS/')

%% 
h = [1 0;0 1];
Nt = 128;

Ey = 2.62e11;
rho = 1280;
Ar = 31.75e-3*54e-3;
% Klib = struct('K', @(w,xi) w/sqrt(Ey/rho), ...
%     'dKdw', @(w,xi) 1/sqrt(Ey/rho), ...
%     'dKdxi', @(w,xi) 0);
Klib = struct('K', @(w,xi) w/sqrt(Ey/rho));
wcomps = [1j 1; -1j 1];

ell = 1;
pcs = [struct('coords', [0;0.28;1/3], 'wcomps', wcomps);
    struct('coords', [1/3;1], 'wcomps', wcomps)];
% bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1], 'dcofsdw', @(w,xi) [0 0], 'dcofsdxi', @(w,xi) [0 0], 'rhs',0);  % The BCs can have w,xi dependence (beam bcs do).
%     struct('i', 5, 'cofs', @(w,xi) [1 1], 'dcofsdw', @(w,xi) [0 0], 'dcofsdxi', @(w,xi) [0 0], 'rhs',0)];
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1]);  % The BCs can have w,xi dependence (beam bcs do).
    struct('i', 5, 'cofs', @(w,xi) [1 1])];
%% Joint
kJ = 1e9;
cJ = 320;
gJ = 1e8;

% Linear Spring-Damper Joint
cofs = @(w,xi) [(1j*Klib.K(w,xi)*Ey*Ar)*[1, -1, 0, 0]+(kJ-1j*cJ*w)*[1 1 -1 -1];1, -1, -1, 1];
dcofsdw = @(w,xi) [(1j*Klib.dKdw(w,xi)*Ey*Ar)*[1, -1, 0, 0]+(-1j*cJ)*[1 1 -1 -1];0, 0, 0, 0];
dcofsdxi = @(w,xi) [(1j*Klib.dKdxi(w,xi)*Ey*Ar)*[1, -1, 0, 0];0, 0, 0, 0];
% linjoints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, 'dcofsdw', dcofsdw, 'dcofsdxi', dcofsdxi);
linjoints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs);

% Fused Joint
fusjoints = struct('type', 2, 'i', 3, 'j', 4);  % Glued Joint

% Nonlinear Joint
cofs = @(w,xi) [1j*[1 -1 0 0];1 -1 -1 1];
% nljoints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
%     'nl', @(Uw) HDUFF(Uw, kJ, cJ, gJ, h, Nt), ...
%     'nlfcofs', @(w,xi) 1./(Klib.K(w,xi)+eps).*[1; 0], ...
%     'dnlfcofsdw', @(w,xi) -Klib.dKdw(w,xi)./(Klib.K(w,xi).^2+eps).*[1; 0], ...
%     'dnlfcofsdxi', @(w,xi) -Klib.dKdxi(w,xi)./(Klib.K(w,xi).^2+eps).*[1; 0], ...
%     'nldcofs', @(w,xi) [1 1 -1 -1]);
nljoints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) QPHDUFF(Uw, kJ, cJ, gJ, h, Nt), ...
    'nlfcofs', @(w,xi) 1./(Klib.K(w,xi)*Ey*Ar).*[1; 0], ...
    'nldcofs', @(w,xi) [1 1 -1 -1]);
%% Excitation
% excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) (1/2/(2j*(Klib.K(w,xi)+eps)*Ey*Ar))*[-1;1], ...
%     'drcofsdw', @(w,xi) (-Klib.dKdw(w,xi)/2/(2j*(Klib.K(w,xi)^2+eps)*Ey*Ar))*[-1;1], ...
%     'drcofsdxi', @(w,xi) (-Klib.dKdxi(w,xi)/2/(2j*(Klib.K(w,xi)^2+eps)*Ey*Ar))*[-1;1]);
excs = struct('i', 2, 'nh', [1 0], 'rcofs', @(w,xi) (1/2/(2j*(Klib.K(w,xi)+eps)*Ey*Ar))*[-1;1]);

%% Preprocess Everything
[~, ~, linjoints, ~, ~] = WBPREPROC(pcs, bcs, linjoints, excs, Klib);
[~, ~, fusjoints, ~, ~] = WBPREPROC(pcs, bcs, fusjoints, excs, Klib);
[pcs, bcs, nljoints, excs, Klib] = WBPREPROC(pcs, bcs, nljoints, excs, Klib);

%% Setup Properties for HB
Famp = 150e5;
Npts = pcs(end).irange(end);
Npcs = length(pcs);
Nwc = 2;
Nh = size(h,1);
Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));

Wst = 78e3;

ari0 = zeros(pcs(end).irange(end)*Nwc*Nhc,1);
ws = Wst*[1;pi];  % Frequencies

%% Conduct HB with full form
opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
aris = fsolve(@(ari) WVHBRESFUN([ari; ws], Famp, h, pcs, bcs, nljoints, Klib), ari0, opt);
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);
acso = zeros(Npts*Nh, 1);
acso([zinds; hinds]) = [aris(rinds0); aris(rinds)+1j*aris(iinds)];

%% Conduct HB with reduced form
ari0r = zeros(length(pcs)*Nwc*Nhc, 1);
opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
arisr = fsolve(@(ari) WVHBRESFUNr([ari; ws], Famp, h, pcs, bcs, nljoints, Klib), ari0r, opt);

[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npcs*Nwc, h);
acsr = zeros(Npcs*Nh, 1);
acsr([zinds; hinds]) = [arisr(rinds0); arisr(rinds)+1j*arisr(iinds)];

[Rh, ~, ~, Ri] = MAPr2COMPS([ws;0],h,pcs,Klib);

acs = Rh*acsr-Ri*Famp;



