clc
clear all

%% 
Ey = 210e9;
rho = 7800;
Ar = 2*pi*1e-2^2;
% Klib = struct('K', @(w,xi) w/sqrt(Ey/rho), ...
%     'dKdw', @(w,xi) 1/sqrt(Ey/rho), ...
%     'dKdxi', @(w,xi) 0);
Klib = struct('K', @(w,xi) w/sqrt(Ey/rho));
wcomps = [1j 1; -1j 1];

ell = 1;
pcs = [struct('coords', [0;ell/2;ell], 'wcomps', wcomps);
    struct('coords', [ell;3*ell/2;2*ell], 'wcomps', wcomps)];
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1]);  % The BCs can have w,xi dependence (beam bcs do).
    struct('i', 6, 'cofs', @(w,xi) [1 1])];
excs = [struct('i', 2, 'nh', [1 0], 'rcofs', @(w,xi) (1/2/(2j*(Klib.K(w,xi)+eps)*Ey*Ar))*[-1;1]);
    struct('i', 5, 'nh', [0 1], 'rcofs', @(w,xi) (1/2/(2j*(Klib.K(w,xi)+eps)*Ey*Ar))*[-1;1])];

%% Joint
% Linear Spring-Damper Joint
kJ = 1e8;
cJ = 400;
cofs = @(w,xi) [(1j*Klib.K(w,xi)*Ey*Ar)*[1, -1, 0, 0]+(kJ-1j*cJ*w)*[1 1 -1 -1];1, -1, -1, 1];
joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs);

% % Fused Joint
% joints = struct('type', 2, 'i', 3, 'j', 4);  % Glued Joint

%% Preprocess
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%% 
h = [1 0;0 1];
Nw = 1000;
Ws = linspace(eps, 1.5e5, Nw);
As = zeros(pcs(end).irange(end)*2*size(h,1), Nw);
for iw = 1:Nw
    Oms = Ws(iw)*[1;pi];
    [Amat, ~, ~, Fv] = WVAMATQP([Oms;0],h,pcs,bcs,joints,Klib);
    As(:, iw) = Amat\Fv;
end

%%
figure(1)
clf()
plot(Ws, abs(2*sum(As((excs(1).i-1)*2+(1:2),:))), '.-'); hold on
plot(Ws*pi, abs(2*sum(As(pcs(end).irange(end)*2+(excs(1).i-1)*2+(1:2),:)))); hold on

set(gca, 'YScale', 'log')