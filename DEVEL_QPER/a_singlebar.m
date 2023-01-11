clc
clear all
addpath('../ROUTINES/WBM/')

%% 
Ey = 210e9;
rho = 7800;
% Klib = struct('K', @(w,xi) w/sqrt(Ey/rho), ...
%     'dKdw', @(w,xi) 1/sqrt(Ey/rho), ...
%     'dKdxi', @(w,xi) 0);
Klib = struct('K', @(w,xi) w/sqrt(Ey/rho));
wcomps = [1j 1; -1j 1];

ell = 1;
Ar = 1e-2;
pcs = struct('coords', [0;ell/3;2*ell/3;ell], 'wcomps', wcomps);
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1], 'rhs',@(w,xi) 0);
    struct('i', 4, 'cofs', @(w,xi) [1 1], 'rhs',@(w,xi) 0)];
joints = [];
excs = [struct('i', 2, 'nh', [0 1], 'rcofs', @(w,xi) (1/2/(2j*(Klib.K(w,xi)+eps)*Ey*Ar))*[-1;1]);
    struct('i', 3, 'nh', [1 0], 'rcofs', @(w,xi) (1/2/(2j*(Klib.K(w,xi)+eps)*Ey*Ar))*[-1;1])];
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%%
h = [1 0;0 1];

Oms = [1e5; pi*1e5];
[Amat, dAmatdw, dAmatdxi, Fv] = WVAMATrQP([Oms;0],h,pcs,bcs,joints,Klib);

%%
Wresa = (1:10)*pi/ell*sqrt(Ey/rho);  % Fix-Fix
% Wresa = (2*(1:10)-1)*pi/2/ell*sqrt(Ey/rho);  % Fix-Free

Nw = 1000;
Ws = linspace(eps, Wresa(end), Nw);
As = zeros(size(Amat,1), Nw);
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