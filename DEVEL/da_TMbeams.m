clc
clear all
addpath('../ROUTINES/SOLVERS/')

%%
Ey = 190e9;
G = 77.5e9;
nu = 0.29;
rho = 7680;
wid = 0.1;
brd = 0.4;
Ar = wid*brd;
Iy = wid^3*brd/12;
L0 = 2.0;
kap = 10*(1+nu)/(12+11*nu);

% Dispersion Relationship
a = Ey*kap*G*Iy;
b = struct('v', @(w, xi) -(w.^2*rho*Ar+rho*Ey*Iy), ...
    'dw', @(w, xi) -(2*w*rho*Ar), ...
    'dxi', @(w, xi) 0);
c = struct('v', @(w, xi) w.^4*rho^2*Ar-w.^2*rho*kap*G*Ar, ...
    'dw', @(w, xi) 4*w.^3*rho^2*Ar-2*w*rho*kap*G*Ar, ...
    'dxi', @(w, xi) 0);
Klib = [struct('K', @(w,xi) sqrt((-b.v(w,xi)+sqrt(b.v(w,xi).^2-4*a*c.v(w,xi)))/(2*a)), ...
    'dKdw', @(w,xi) 0.5*((-b.v(w,xi)+sqrt(b.v(w,xi).^2-4*a*c.v(w,xi)))/(2*a)).^(-0.5).*(-b.dw(w,xi)+0.5*(b.v(w,xi).^2-4*a*c.v(w,xi)).^(-0.5).*(2*b.v(w,xi).*b.dw(w,xi)-4*a*c.dw(w,xi)))/(2*a), ...
    'dKdxi', @(w,xi) 0.5*((-b.v(w,xi)+sqrt(b.v(w,xi).^2-4*a*c.v(w,xi)))/(2*a)).^(-0.5).*(-b.dxi(w,xi)+0.5*(b.v(w,xi).^2-4*a*c.v(w,xi)).^(-0.5).*(2*b.v(w,xi).*b.dxi(w,xi)-4*a*c.dxi(w,xi)))/(2*a)); ...
    struct('K', @(w,xi) sqrt((-b.v(w,xi)-sqrt(b.v(w,xi).^2-4*a*c.v(w,xi)))/(2*a)), ...
    'dKdw', @(w,xi) 0.5*((-b.v(w,xi)-sqrt(b.v(w,xi).^2-4*a*c.v(w,xi)))/(2*a)).^(-0.5).*(-b.dw(w,xi)-0.5*(b.v(w,xi).^2-4*a*c.v(w,xi)).^(-0.5).*(2*b.v(w,xi).*b.dw(w,xi)-4*a*c.dw(w,xi)))/(2*a), ...
    'dKdxi', @(w,xi) 0.5*((-b.v(w,xi)-sqrt(b.v(w,xi).^2-4*a*c.v(w,xi)))/(2*a)).^(-0.5).*(-b.dxi(w,xi)-0.5*(b.v(w,xi).^2-4*a*c.v(w,xi)).^(-0.5).*(2*b.v(w,xi).*b.dxi(w,xi)-4*a*c.dxi(w,xi)))/(2*a))];
wcomps = [1j 1;  % +jK1
         -1j 1;  % -jK1
          1j 2;  % +jK2
         -1j 2]; % -jK2

% The following represents the ratio of "theta"'s wave coefficient wrt
% "v"'s, i.e., the ratio of rotation wave coefficient w.r.t. translational
% wave coefficient
rel = struct('P', @(w,xi,cc,K) ((cc*K.K(w,xi)).^2*kap*G*Ar-w.^2*rho*Ar)./(-1j*cc*K.K(w,xi)*kap*G*Ar), ...
    'dPdw', @(w,xi,cc,K) ((-1j*cc*K.K(w,xi)*kap*G*Ar).*(2*cc*K.K(w,xi).*K.dKdw(w,xi)*kap*G*Ar-2*w*rho*Ar)-((cc*K.K(w,xi)).^2*kap*G*Ar-w.^2*rho*Ar).*(-1j*cc*K.dKdw(w,xi)*kap*G*Ar))./(-1j*cc*K.K(w,xi)*kap*G*Ar).^2, ...
    'dPdxi', @(w,xi,cc,K) ((-1j*cc*K.K(w,xi)*kap*G*Ar).*(2*cc*K.K(w,xi).*K.dKdxi(w,xi)*kap*G*Ar)-((cc*K.K(w,xi)).^2*kap*G*Ar-w.^2*rho*Ar).*(-1j*cc*K.dKdxi(w,xi)*kap*G*Ar))./(-1j*cc*K.K(w,xi)*kap*G*Ar).^2);

% Pieces
pcs = struct('coords', [0; L0]);
% bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1], 'dcofsdw', @(w,xi) [0 0 0 0], 'dcofsdxi', @(w,xi) [0 0 0 0]);
%     struct('i', 1, 'cofs', @(w,xi) [rel.P(w,xi,1,Klib(1)) rel.P(w,xi,-1,Klib(1)) rel.P(w,xi,1,Klib(2)) rel.P(w,xi,-1,Klib(2))], ...
%     'dcofsdw', @(w,xi) [rel.dPdw(w,xi,1,Klib(1)) rel.dPdw(w,xi,-1,Klib(1)) rel.dPdw(w,xi,1,Klib(2)) rel.dPdw(w,xi,-1,Klib(2))], ...
%     'dcofsdxi', @(w,xi) [rel.dPdxi(w,xi,1,Klib(1)) rel.dPdxi(w,xi,-1,Klib(1)) rel.dPdxi(w,xi,1,Klib(2)) rel.dPdxi(w,xi,-1,Klib(2))]);
%     struct('i', 2, 'cofs', @(w,xi) [1 1 1 1], 'dcofsdw', @(w,xi) [0 0 0 0], 'dcofsdxi', @(w,xi) [0 0 0 0]);
%     struct('i', 2, 'cofs', @(w,xi) [rel.P(w,xi,1,Klib(1)) rel.P(w,xi,-1,Klib(1)) rel.P(w,xi,1,Klib(2)) rel.P(w,xi,-1,Klib(2))], ...
%     'dcofsdw', @(w,xi) [rel.dPdw(w,xi,1,Klib(1)) rel.dPdw(w,xi,-1,Klib(1)) rel.dPdw(w,xi,1,Klib(2)) rel.dPdw(w,xi,-1,Klib(2))], ...
%     'dcofsdxi', @(w,xi) [rel.dPdxi(w,xi,1,Klib(1)) rel.dPdxi(w,xi,-1,Klib(1)) rel.dPdxi(w,xi,1,Klib(2)) rel.dPdxi(w,xi,-1,Klib(2))])];
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1], 'dcofsdw', @(w,xi) [0 0 0 0], 'dcofsdxi', @(w,xi) [0 0 0 0]);
    struct('i', 1, 'cofs', @(w,xi) [rel.P(w,xi,1,Klib(1)) rel.P(w,xi,-1,Klib(1)) rel.P(w,xi,1,Klib(2)) rel.P(w,xi,-1,Klib(2))], ...
    'dcofsdw', @(w,xi) [rel.dPdw(w,xi,1,Klib(1)) rel.dPdw(w,xi,-1,Klib(1)) rel.dPdw(w,xi,1,Klib(2)) rel.dPdw(w,xi,-1,Klib(2))], ...
    'dcofsdxi', @(w,xi) [rel.dPdxi(w,xi,1,Klib(1)) rel.dPdxi(w,xi,-1,Klib(1)) rel.dPdxi(w,xi,1,Klib(2)) rel.dPdxi(w,xi,-1,Klib(2))]);
    struct('i', 2, 'cofs', @(w,xi) [1j*Klib(1).K(w,xi)-rel.P(w,xi,1,Klib(1)) -1j*Klib(1).K(w,xi)-rel.P(w,xi,-1,Klib(1)) 1j*Klib(2).K(w,xi)-rel.P(w,xi,1,Klib(2)) -1j*Klib(2).K(w,xi)-rel.P(w,xi,-1,Klib(2))], ...
    'dcofsdw', @(w,xi) [1j*Klib(1).dKdw(w,xi)-rel.dPdw(w,xi,1,Klib(1)) -1j*Klib(1).dKdw(w,xi)-rel.dPdw(w,xi,-1,Klib(1)) 1j*Klib(2).dKdw(w,xi)-rel.dPdw(w,xi,1,Klib(2)) -1j*Klib(2).dKdw(w,xi)-rel.dPdw(w,xi,-1,Klib(2))], ...
    'dcofsdxi', @(w,xi) [1j*Klib(1).dKdxi(w,xi)-rel.dPdxi(w,xi,1,Klib(1)) -1j*Klib(1).dKdxi(w,xi)-rel.dPdxi(w,xi,-1,Klib(1)) 1j*Klib(2).dKdxi(w,xi)-rel.dPdxi(w,xi,1,Klib(2)) -1j*Klib(2).dKdxi(w,xi)-rel.dPdxi(w,xi,-1,Klib(2))]);
    struct('i', 2, 'cofs', @(w,xi) [Klib(1).K(w,xi)*rel.P(w,xi,1,Klib(1)) -Klib(1).K(w,xi)*rel.P(w,xi,-1,Klib(1)) Klib(2).K(w,xi)*rel.P(w,xi,1,Klib(2)) -Klib(2).K(w,xi)*rel.P(w,xi,-1,Klib(2))], ...
    'dcofsdw', @(w,xi) [Klib(1).dKdw(w,xi)*rel.P(w,xi,1,Klib(1))+Klib(1).K(w,xi)*rel.dPdw(w,xi,1,Klib(1)) -Klib(1).dKdw(w,xi)*rel.P(w,xi,-1,Klib(1))-Klib(1).K(w,xi)*rel.dPdw(w,xi,1,Klib(1)) Klib(2).dKdw(w,xi)*rel.P(w,xi,1,Klib(2))+Klib(2).K(w,xi)*rel.dPdw(w,xi,1,Klib(2)) -Klib(2).dKdw(w,xi)*rel.P(w,xi,-1,Klib(2))-Klib(2).K(w,xi)*rel.dPdw(w,xi,-1,Klib(2))], ...
    'dcofsdxi', @(w,xi) [Klib(1).dKdxi(w,xi)*rel.P(w,xi,1,Klib(1))+Klib(1).K(w,xi)*rel.dPdxi(w,xi,1,Klib(1)) -Klib(1).dKdxi(w,xi)*rel.P(w,xi,-1,Klib(1))-Klib(1).K(w,xi)*rel.dPdxi(w,xi,1,Klib(1)) Klib(2).dKdxi(w,xi)*rel.P(w,xi,1,Klib(2))+Klib(2).K(w,xi)*rel.dPdxi(w,xi,1,Klib(2)) -Klib(2).dKdxi(w,xi)*rel.P(w,xi,-1,Klib(2))-Klib(2).K(w,xi)*rel.dPdxi(w,xi,-1,Klib(2))])];

% Joints
joints = [];

% Excitation
excs = [];

%% Preprocessing
[pcs, bcs, joints, excs] = WBPREPROC(pcs, bcs, joints, excs, wcomps);

%% Amatrix
iw = 1;
h = 1;

Nw = 2000;
Ws = linspace(eps, 1e4, Nw);
Dv = zeros(size(Ws));
for iw=1:Nw
    Amat = WVAMAT([Ws(iw);0],h,pcs,bcs,joints,Klib,wcomps);
    Dv(iw) = det([real([Amat 1j*Amat]); imag([Amat 1j*Amat])]);
end

% %% Analytical solution
% Wres = 2*pi*3.56*sqrt(Ey*Iy/(rho*Ar*L0^4));
% 
% fan = @(w) cos(sqrt(w)*sqrt(sqrt(rho*Ar/Ey/Iy))*L0).*cosh(sqrt(w)*sqrt(sqrt(rho*Ar/Ey/Iy))*L0)-1;
% fz = @(lam) deal(cos(lam*L0)*cosh(lam*L0)-1, -L0*sin(lam*L0)*cosh(lam*L0)+cos(lam*L0)*sinh(lam*L0));
% opt = optimoptions('fsolve', 'SpecifyObjectiveGradient',true, 'Display','iter');
% Lamsol = fsolve(fz, sqrt(Wres)*(rho*Ar/(Ey*Iy))^(0.25), opt);

%%
ebKlib = struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25), ...
    'dKdw', @(w,xi) 0.5/sqrt(w)*(rho*Ar/Ey/Iy)^(0.25), ...
    'dKdxi', @(w,xi) 0);
fixfix = @(w) 2*(1-cos(ebKlib.K(w,0)*L0).*cosh(ebKlib.K(w,0)*L0));
fixfree = @(w) 2*(1+cos(ebKlib.K(w,0)*L0).*cosh(ebKlib.K(w,0)*L0));

figure(1)
clf()
semilogy(Ws, Dv); hold on
% plot(Ws, abs(fixfix(Ws)))
plot(Ws, abs(fixfree(Ws)))