clc
clear all
addpath('../ROUTINES/SOLVERS/')

%% 
h = 1:15;
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
    'nl', @(Uw) HDUFF(Uw, kJ, cJ, gJ, h, Nt), ...
    'nlfcofs', @(w,xi) 1./(Klib.K(w,xi)*Ey*Ar).*[1; 0], ...
    'nldcofs', @(w,xi) [1 1 -1 -1]);
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
Npcs = length(pcs);
Nwc = 2;
Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));

Wst = 78e3;
Wen = 100e3;
dw = 0.1;

Copt = struct('Nmax', 100, 'angopt', 2e-1, 'DynDscale', 1);

%% HB Fresp - Reduced form
Hs = {1, 1:3, 1:6, 1:12, 1:24, 1:48, 1:96};
Nhs = cellfun(@(h) length(h), Hs);
Ttks1 = zeros(size(Hs));
% tic
for hi=1:length(Hs)
    nljoints.nl = @(Uw) HDUFF(Uw, kJ, cJ, gJ, Hs{hi}, Nt);
    Nhc = sum((Hs{hi}==0)+2*(Hs{hi}~=0));
    ari0 = zeros(Npcs*Nwc*Nhc,1);
    tic
    ariwC = CONTINUE(@(ariw) WVHBRESFUNr(ariw, Famp, Hs{hi}, pcs, bcs, nljoints, Klib),  ...
        ari0, Wst, Wen, dw, Copt);
    Ttks1(hi) = toc;

    fprintf('Done %d\n', hi);
end
% tk1=toc;

[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npcs*Nwc, h);
acC1 = zeros(Npcs*Nwc*Nh+1, size(ariwC,2));
acC1([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];

%% HB Fresp - Full form
Ttks2 = zeros(size(Hs));

for hi=1:length(Hs)
    nljoints.nl = @(Uw) HDUFF(Uw, kJ, cJ, gJ, Hs{hi}, Nt);
    Nhc = sum((Hs{hi}==0)+2*(Hs{hi}~=0));
    ari0 = zeros(Npts*Nwc*Nhc,1);

    tic
    ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famp, Hs{hi}, pcs, bcs, nljoints, Klib),  ...
        zeros(Npts*Nwc*Nhc,1), Wst, Wen, dw, Copt);
    Ttks2(hi)=toc;

    fprintf('Done %d\n', hi);
end

[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);
acC2 = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
acC2([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];
%%
figure(1)
clf()
loglog(Nhs, [Ttks1;Ttks2], 'o-', 'LineWidth', 2); 
grid on
xlabel('Number of Harmonics')
ylabel('Computation Time')
legend('Naive Approach', 'Eliminated Approach')


opi1 = 3:4;
opi2 = 9:10;
figure(2)
clf()
plot(acC1(end,:)/1e3, abs(2*sum(acC1(opi1,:))), '.-'); hold on
plot(acC2(end,:)/1e3, abs(2*sum(acC2(opi2,:))), 'o'); hold on
xlim([80 81]);
ylim([0.05 0.3])
xlabel('Frequency (k rad/s)')
ylabel('Response (m)')