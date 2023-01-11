clc
clear all
addpath('../ROUTINES/SOLVERS/')

%%
Ey = 190e9;
G = 77.5e9;
nu = 0.29;
rho = 7680;
wid = 0.2;
brd = 0.4;
Ar = wid*brd;
Iy = wid^3*brd/12;
L0 = 2.0;
kap = 10*(1+nu)/(12+11*nu);
% Klib = struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25), ...
%     'dKdw', @(w,xi) 0.5/sqrt(w)*(rho*Ar/Ey/Iy)^(0.25), ...
%     'dKdxi', @(w,xi) 0);
Klib = struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25));
wcomps = [1 1;
         -1 1;
         1j 1;
        -1j 1];

pcs = struct('coords', [0;L0], 'wcomps', wcomps);
% % fix-fix
% bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1], 'dcofsdw', [], 'dcofsdxi', []);
%     struct('i', 1, 'cofs', @(w,xi) [1 -1 1j -1j], 'dcofsdw', [], 'dcofsdxi', []);
%     struct('i', 2, 'cofs', @(w,xi) [1 1 1 1], 'dcofsdw', [], 'dcofsdxi', []);
%     struct('i', 2, 'cofs', @(w,xi) [1 -1 1j -1j], 'dcofsdw', [], 'dcofsdxi', [])];
% % fix-free
% bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1], 'dcofsdw', [], 'dcofsdxi', []);
%     struct('i', 1, 'cofs', @(w,xi) [1 -1 1j -1j], 'dcofsdw', [], 'dcofsdxi', []);
%     struct('i', 2, 'cofs', @(w,xi) [1 1 -1 -1], 'dcofsdw', [], 'dcofsdxi', []);
%     struct('i', 2, 'cofs', @(w,xi) [1 -1 -1j 1j], 'dcofsdw', [], 'dcofsdxi', [])];
% fix-fix
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1]);
    struct('i', 1, 'cofs', @(w,xi) [1 -1 1j -1j]);
    struct('i', 2, 'cofs', @(w,xi) [1 1 1 1]);
    struct('i', 2, 'cofs', @(w,xi) [1 -1 1j -1j])];
% fix-free
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1]);
    struct('i', 1, 'cofs', @(w,xi) [1 -1 1j -1j]);
    struct('i', 2, 'cofs', @(w,xi) [1 1 -1 -1]);
    struct('i', 2, 'cofs', @(w,xi) [1 -1 -1j 1j])];

joints = [];
excs = [];

[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%% Amatrix
iw = 1;
h = 1;

Nw = 1000;
Ws = linspace(eps, 1e4, Nw);
Dv = zeros(size(Ws));
for iw=1:Nw
    Amat = WVAMAT([Ws(iw);0],h,pcs,bcs,joints,Klib);
    Dv(iw) = det([real([Amat 1j*Amat]); imag([Amat 1j*Amat])]);
end

%% Analytical
fixfix = @(w) 2*(1-cos(Klib.K(w,0)*L0).*cosh(Klib.K(w,0)*L0));
fixfree = @(w) 2*(1+cos(Klib.K(w,0)*L0).*cosh(Klib.K(w,0)*L0));

figure(1)
clf()
% semilogy(Ws, [Dv; abs(fixfix(Ws))])
semilogy(Ws, [Dv; abs(fixfree(Ws))])