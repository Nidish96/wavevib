clc
clear all
addpath('../ROUTINES/SOLVERS/')

%% 
h = [1];
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
    struct('coords', [1/3;1], 'wcomps', wcomps);
    struct('coords', [1;4/3], 'wcomps', wcomps)];
% bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1], 'dcofsdw', @(w,xi) [0 0], 'dcofsdxi', @(w,xi) [0 0], 'rhs',0);  % The BCs can have w,xi dependence (beam bcs do).
%     struct('i', 7, 'cofs', @(w,xi) [1 1], 'dcofsdw', @(w,xi) [0 0], 'dcofsdxi', @(w,xi) [0 0], 'rhs',0)];
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1]);  % The BCs can have w,xi dependence (beam bcs do).
    struct('i', 7, 'cofs', @(w,xi) [1 1])];
%% Joint
kJ = 1e9;
cJ = 320;
gJ = 1e8;

% Linear Spring-Damper Joint
cofs = @(w,xi) [(1j*Klib.K(w,xi)*Ey*Ar)*[1, -1, 0, 0]+(kJ-1j*cJ*w)*[1 1 -1 -1];1, -1, -1, 1];
dcofsdw = @(w,xi) [(1j*Klib.dKdw(w,xi)*Ey*Ar)*[1, -1, 0, 0]+(-1j*cJ)*[1 1 -1 -1];0, 0, 0, 0];
dcofsdxi = @(w,xi) [(1j*Klib.dKdxi(w,xi)*Ey*Ar)*[1, -1, 0, 0];0, 0, 0, 0];
% linjoints = [struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, 'dcofsdw', dcofsdw, 'dcofsdxi', dcofsdxi);
%     struct('type', 2, 'i', 5, 'j', 6, 'cofs', cofs, 'dcofsdw', dcofsdw, 'dcofsdxi', dcofsdxi)];
linjoints = [struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs);
    struct('type', 2, 'i', 5, 'j', 6, 'cofs', cofs)];

% Fused Joint
fusjoints = [struct('type', 2, 'i', 3, 'j', 4);
    struct('type', 2, 'i', 5, 'j', 6)];  % Glued Joint

% Nonlinear Joint
cofs = @(w,xi) [1j*Ey*Ar*[1 -1 0 0];1 -1 -1 1];
% nljoints = [struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
%     'nl', @(Uw) HDUFF(Uw, kJ, cJ, gJ, h, Nt), ...
%     'nlfcofs', @(w,xi) 1./(Klib.K(w,xi)+eps).*[1; 0], ...
%     'dnlfcofsdw', @(w,xi) -Klib.dKdw(w,xi)./(Klib.K(w,xi).^2+eps).*[1; 0], ...
%     'dnlfcofsdxi', @(w,xi) -Klib.dKdxi(w,xi)./(Klib.K(w,xi).^2+eps).*[1; 0], ...
%     'nldcofs', @(w,xi) [1 1 -1 -1]);
%     struct('type', 2, 'i', 5, 'j', 6, 'cofs', cofs, ...
%     'nl', @(Uw) HDUFF(Uw, kJ, cJ, gJ, h, Nt), ...
%     'nlfcofs', @(w,xi) 1./(Klib.K(w,xi)+eps).*[1; 0], ...
%     'dnlfcofsdw', @(w,xi) -Klib.dKdw(w,xi)./(Klib.K(w,xi).^2+eps).*[1; 0], ...
%     'dnlfcofsdxi', @(w,xi) -Klib.dKdxi(w,xi)./(Klib.K(w,xi).^2+eps).*[1; 0], ...
%     'nldcofs', @(w,xi) [1 1 -1 -1])];
nljoints = [struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJ, cJ, gJ, h, Nt), ...
    'nlfcofs', @(w,xi) 1./(Klib.K(w,xi)+eps).*[1; 0], ...
    'nldcofs', @(w,xi) [1 1 -1 -1]);
    struct('type', 2, 'i', 5, 'j', 6, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJ, cJ, gJ, h, Nt), ...
    'nlfcofs', @(w,xi) 1./(Klib.K(w,xi)+eps).*[1; 0], ...
    'nldcofs', @(w,xi) [1 1 -1 -1])];
%% Excitation
% excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) (1/2/(2j*(Klib.K(w,xi)+eps)*Ey*Ar))*[-1;1], ...
%     'drcofsdw', @(w,xi) (-Klib.dKdw(w,xi)/2/(2j*(Klib.K(w,xi)^2+eps)*Ey*Ar))*[-1;1], ...
%     'drcofsdxi', @(w,xi) (-Klib.dKdxi(w,xi)/2/(2j*(Klib.K(w,xi)^2+eps)*Ey*Ar))*[-1;1]);
excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) (1/2/(2j*(Klib.K(w,xi)+eps)*Ey*Ar))*[-1;1]);

%% Preprocess Everything
[~, ~, linjoints, ~, ~] = WBPREPROC(pcs, bcs, linjoints, excs, Klib);
[~, ~, fusjoints, ~, ~] = WBPREPROC(pcs, bcs, fusjoints, excs, Klib);
[pcs, bcs, nljoints, excs, Klib] = WBPREPROC(pcs, bcs, nljoints, excs, Klib);

%% Conduct HB
Famp = 150e5;
Npts = pcs(end).irange(end);
Nwc = 2;
Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));

Wst = 25e3;
Wen = 26e3;
dw = 0.1;

[Amat, ~, ~, Fv] = WVAMAT([Wst;0],h,pcs,bcs,linjoints,Klib,'r');
ari0 = Amat\Fv*Famp;

Copt = struct('Nmax', 100, 'angopt', 2e-1, 'DynDscale', 1);
ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famp, h, pcs, bcs, nljoints, Klib),  ...
    ari0, Wst, Wen, dw, Copt);

[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);
acC = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
acC([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];

%% Redo with original complex mx
aciC = zeros(size(acC));
for iw=1:size(acC,2)
    [Amat, ~, ~, Fv] = WVAMAT([acC(end,iw);0],h,pcs,bcs,linjoints,Klib);
    aciC(1:end-1,iw) = Amat\Fv*Famp;
    aciC(end,iw) = acC(end,iw);
end
%%
opi =  9:10;
figure(2)
% clf()
% semilogy(Ws/1e3, 2*abs(sum(Us(opi,:)))/Famp, '.-'); hold on
% semilogy(Ws/1e3, 2*abs(sum(Unlsc(opi,:)))/Famp, '.-'); hold on
plot(acC(end,:)/1e3, 2*abs(sum(acC(opi,:))), '.-'); hold on
% xlim([80 81]);
% ylim([0.05 0.3])