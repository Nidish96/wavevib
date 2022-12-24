clc
clear all

%% 
Ey = 210e9;
rho = 7800;
% Klib = struct('K', @(w,xi) w/sqrt(Ey/rho), ...
%     'dKdw', @(w,xi) 1/sqrt(Ey/rho), ...
%     'dKdxi', @(w,xi) 0);
Klib = struct('K', @(w,xi) w/sqrt(Ey/rho));
wcomps = [1j 1; -1j 1];

ell = 1;
pcs = struct('coords', [0;ell/2;ell], 'wcomps', wcomps);
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1], 'rhs',@(w,xi) 0);
    struct('i', 3, 'cofs', @(w,xi) [1 1], 'rhs',@(w,xi) 0)];
joints = [];
excs = [];
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%%
h = [1];
[Amat, dAmatdw, dAmatdxi] = WVAMAT([1e5;0],h,pcs,bcs,joints,Klib);

%%
Wresa = (1:10)*pi/ell*sqrt(Ey/rho);  % Fix-Fix
% Wresa = (2*(1:10)-1)*pi/2/ell*sqrt(Ey/rho);  % Fix-Free

Nw = 1000;
Ws = linspace(0, Wresa(end), Nw);
Dv = zeros(1,Nw);
for iw = 1:Nw
    Amat = WVAMAT([Ws(iw);0],h,pcs,bcs,joints,Klib);
    Dv(iw) = det([real([Amat 1j*Amat]); imag([Amat 1j*Amat])]);
end

%%
figure(1)
clf()
plot(Ws, Dv); hold on
for ri=1:10
    plot(Wresa(ri)*[1 1], ylim, 'k-')
end