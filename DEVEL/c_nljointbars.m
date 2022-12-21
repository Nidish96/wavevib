clc
clear all

%% 
h = [1];
Nt = 128;

Ey = 2.62e11;
rho = 1280;
Ar = 31.75e-3*54e-3;
Klib = struct('K', @(w,xi) w/sqrt(Ey/rho), ...
    'dKdw', @(w,xi) 1/sqrt(Ey/rho), ...
    'dKdxi', @(w,xi) 0);
wcomps = [1j 1; -1j 1];

ell = 1;
pcs = [struct('coords', [0;0.28;1/3]);
    struct('coords', [1/3;1])];
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1], 'dcofsdw', @(w,xi) [0 0], 'dcofsdxi', @(w,xi) [0 0], 'rhs',0);  % The BCs can have w,xi dependence (beam bcs do).
    struct('i', 5, 'cofs', @(w,xi) [1 1], 'dcofsdw', @(w,xi) [0 0], 'dcofsdxi', @(w,xi) [0 0], 'rhs',0)];
%% Joint
kJ = 1e9;
cJ = 320;
gJ = 1e8;

% Linear Spring-Damper Joint
cofs = @(w,xi) [(1j*Klib.K(w,xi)*Ey*Ar)*[1, -1, 0, 0]+(kJ-1j*cJ*w)*[1 1 -1 -1];1, -1, -1, 1];
dcofsdw = @(w,xi) [(1j*Klib.dKdw(w,xi)*Ey*Ar)*[1, -1, 0, 0]+(-1j*cJ)*[1 1 -1 -1];0, 0, 0, 0];
dcofsdxi = @(w,xi) [(1j*Klib.dKdxi(w,xi)*Ey*Ar)*[1, -1, 0, 0];0, 0, 0, 0];
linjoints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, 'dcofsdw', dcofsdw, 'dcofsdxi', dcofsdxi);

% Fused Joint
fusjoints = struct('type', 2, 'i', 3, 'j', 4);  % Glued Joint

% Nonlinear Joint
cofs = @(w,xi) [1 -1 0 0;1 -1 -1 1];
nljoints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJ, cJ, gJ, h, Nt), ...
    'nlfcofs', @(w,xi) 1/(1j*Klib.K(w,xi)*Ey*Ar)*[1; 0], ...
    'dnlfcofsdw', @(w,xi) -Klib.dKdw(w,xi)/(1j*Klib.K(w,xi)^2*Ey*Ar)*[1; 0], ...
    'dnlfcofsdxi', @(w,xi) -Klib.dKdxi(w,xi)/(1j*Klib.K(w,xi)^2*Ey*Ar)*[1; 0], ...
    'nldcofs', @(w,xi) [1 1 -1 -1]);

joints = nljoints;
% joints = linjoints;
%% Excitation
excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) (1/2/(2j*Klib.K(w,xi)*Ey*Ar))*[-1;1], ...
    'drcofsdw', @(w,xi) (-Klib.dKdw(w,xi)/2/(2j*Klib.K(w,xi)^2*Ey*Ar))*[-1;1], ...
    'drcofsdxi', @(w,xi) (-Klib.dKdxi(w,xi)/2/(2j*Klib.K(w,xi)^2*Ey*Ar))*[-1;1]);

%% Preprocess Everything
[pcs, bcs, joints, excs] = WBPREPROC(pcs, bcs, joints, excs, wcomps);

%% Conduct HB
Famp = 150e5*2;
Npts = pcs(end).irange(end);
Nwc = 2;
Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
Nw = 500;
% Ws = linspace(1.4e4, 2e4, Nw);
Ws = linspace(79e3, 82e3, Nw);

opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'off');
Unls = zeros(pcs(end).irange(end)*2*Nhc, Nw);
ari0 = zeros(Npts*Nwc*Nhc,1);
for iw=1:Nw
    Unls(:,iw) = fsolve(@(ari) WVHBRESFUN([ari; Ws(iw)], Famp, h, pcs, bcs, joints, Klib, wcomps), ari0, opt);
    ari0 = Unls(:,iw);
    fprintf('%d/%d\n', iw,Nw);
end
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);
Unlsc = zeros(Npts*Nwc*Nh, Nw);
Unlsc([zinds hinds], :) = [Unls(rinds0,:); Unls(rinds,:)+1j*Unls(iinds,:)];

%%
opi =  9:10;
figure(2)
% clf()
% semilogy(Ws, 2*abs(sum(Us(opi,:)))/Famp, '.-'); hold on
semilogy(Ws, 2*abs(sum(Unlsc(opi,:)))/Famp, '.-'); hold on
xl = xlim;
xlim(xl);