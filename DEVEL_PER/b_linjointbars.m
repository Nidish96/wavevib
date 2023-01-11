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
%% Joint
% Linear Spring-Damper Joint
kJ = 1e8;
cJ = 400;
% cofs = @(w,xi) [(1j*Klib.K(w,xi)*Ey*Ar)*[1, -1, 0, 0]+(kJ-1j*cJ*w)*[1 1 -1 -1];1, -1, -1, 1];
% dcofsdw = @(w,xi) [(1j*Klib.dKdw(w,xi)*Ey*Ar)*[1, -1, 0, 0]+(-1j*cJ)*[1 1 -1 -1];0, 0, 0, 0];
% dcofsdxi = @(w,xi) [(1j*Klib.dKdxi(w,xi)*Ey*Ar)*[1, -1, 0, 0];0, 0, 0, 0];
% joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, 'dcofsdw', dcofsdw, 'dcofsdxi', dcofsdxi);
cofs = @(w,xi) [(1j*Klib.K(w,xi)*Ey*Ar)*[1, -1, 0, 0]+(kJ-1j*cJ*w)*[1 1 -1 -1];1, -1, -1, 1];
joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs);

% % Fused Joint
joints = struct('type', 2, 'i', 3, 'j', 4);  % Glued Joint

% 'type' : 
% 
% 1(Single point being connected to ground at location). 'i' is before and
% 'j' is after joint.
%
% 2(2pts being connected, with relative displacement governing joint), etc.
%% Excitation
% excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) (1/2/(2j*Klib.K(w,xi)*Ey*Ar))*[-1;1], ...
%     'drcofsdw', @(w,xi) (-Klib.dKdw(w,xi)/2/(2j*Klib.K(w,xi)^2*Ey*Ar))*[-1;1], ...
%     'drcofsdxi', @(w,xi) (-Klib.dKdxi(w,xi)/2/(2j*Klib.K(w,xi)^2*Ey*Ar))*[-1;1]);
excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) (1/2/(2j*Klib.K(w,xi)*Ey*Ar))*[-1;1]);

%% Preprocess Everything
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%%
h = [1];
[Amat, dAmatdw, dAmatdxi, Fv, dFvdw, dFvdxi, JEV] = WVAMAT([1e5;0],h,pcs,bcs,joints,Klib);

%%
Wresa = (1:10)*pi/(2*ell)*sqrt(Ey/rho);  % Fix-Fix
% Wresa = (2*(1:10)-1)*pi/2/ell*sqrt(Ey/rho);  % Fix-Free

Nw = 1000;
Ws = linspace(1e-5, Wresa(end), Nw);
Dv = zeros(1,Nw);
Us = zeros(pcs(end).irange(end)*2*length(h), Nw);
for iw = 1:Nw
%     [Amat, ~, ~, Fv] = WVAMAT([Ws(iw);0],h,pcs,bcs,joints,Klib);
%     Dv(iw) = det([real([Amat 1j*Amat]); imag([Amat 1j*Amat])]);
    Dv(iw) = WVLDETFUN([Ws(iw);0],h,pcs,bcs,joints,Klib);
    Us(:, iw) = Amat\Fv;
end

%%
figure(1)
clf()
plot(Ws, abs(Dv)); hold on
for ri=1:10
    plot(Wresa(ri)*[1 1], ylim, 'k-')
end

%%
opi =  3:4;
figure(2)
clf()
semilogy(Ws, 2*abs(sum(Us(opi,:))), '.-'); hold on
for ri=1:10
    plot(Wresa(ri)*[1 1], ylim, 'k-')
end

%%
[Amat, dAmatdw, dAmatdxi, Fv] = WVAMAT([1e5;0],h,pcs,bcs,joints,Klib);