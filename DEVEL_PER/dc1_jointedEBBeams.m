clc
clear all
addpath('../ROUTINES/SOLVERS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

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
Klib = struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25));
wcomps = [1 1;
         -1 1;
         1j 1;
        -1j 1];

pcs = [struct('coords', [0;L0/3;L0], 'wcomps', wcomps);
    struct('coords', [L0;2*L0], 'wcomps', wcomps)];
% fix-fix
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1]);
    struct('i', 1, 'cofs', @(w,xi) [1 -1 1j -1j]);
    struct('i', 5, 'cofs', @(w,xi) [1 1 1 1]);
    struct('i', 5, 'cofs', @(w,xi) [1 -1 1j -1j])];
% % fix-free
% bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1]);
%     struct('i', 1, 'cofs', @(w,xi) [1 -1 1j -1j]);
%     struct('i', 2, 'cofs', @(w,xi) [1 1 -1 -1]);
%     struct('i', 2, 'cofs', @(w,xi) [1 -1 -1j 1j])];

%% Excitation
Mx = @(w,xi) inv([Ey*Iy*Klib.K(w,xi)^3*[-1 1 1j -1j];
    Ey*Iy*Klib.K(w,xi)^2*[-1 -1 1 1];
    [1 1 1 1];
    Klib.K(w,xi)*[1 -1 1j -1j]]);

excs = struct('i', 2, 'nh', 1, ...
    'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0]);

%% Joints
% joints = struct('type', 2, 'i', 3, 'j', 4);

% Nonlinear Joint
h = [1:3];
Nt = 128;
kJs = diag([1e9 1e9]);
cJs = diag([320 320]);
gJs = diag([1e8 0]);
cofs = @(w,xi) [-Klib.K(w,xi)^3*[1 -1 -1j 1j 0 0 0 0];
    -Klib.K(w,xi)^2*[1 1 -1 -1 0 0 0 0];
    -Klib.K(w,xi)^3*[1 -1 -1j 1j -1 1 1j -1j];
    -Klib.K(w,xi)^2*[1 1 -1 -1 -1 -1 1 1]];
joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJs, cJs, gJs, h, Nt), ...
    'nldcofs', @(w,xi) [1 1 1 1 -1 -1 -1 -1; Klib.K(w,xi)*[1 -1 1j -1j -1 1 -1j 1j]], ...
    'nlfcofs', @(w,xi) [eye(2);zeros(2)]/(Ey*Iy));

%% Pre Processing
tic
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);
toc

%% HB
Famp = 1e5; %[1, 10, 20]

Npcs = length(pcs);
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));

Wst = 1055;
Wen = 1075;
dw = 1;

Copt = struct('Nmax', 1000, 'angopt', 1e-1, 'DynDscale', 1);
% fmuls = [1 10 20];
% acC = cell(size(fmuls));

%%
    [Amat, ~, ~, Fv] = WVAMATr([Wst;0], h, pcs, bcs, joints, Klib, 'r');
    ari0 = Amat\Fv*Famp; 
    Copt.Dscale = [abs(ari0+1e-8);Wst];
    
    tic
    ariwC = CONTINUE(@(ariw) WVHBRESFUNr(ariw, Famp, h, pcs, bcs, joints, Klib), ...
        ari0, Wst, Wen, dw, Copt);
    ttk1=toc;
    ttk1 = ttk1/size(ariwC,2);
    [zinds,hinds,rinds0,rinds,iinds] = HINDS(Npcs*Nwc, h);
    acC1 = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
    acC1([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];

% %%
%     [Amat, ~, ~, Fv] = WVAMAT([Wst;0], h, pcs, bcs, joints, Klib, 'r');
%     ari0 = Amat\Fv*Famp;
%     Copt.Dscale = [abs(ari0+1e-6);Wst];
%     
%     tic
%     ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famp, h, pcs, bcs, joints, Klib), ...
%         ari0, Wst, Wen, dw, Copt);
%     ttk2=toc;
%     ttk2 = ttk2/size(ariwC,2);
%     [zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);
%     acC2 = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
%     acC2([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];

%%
opi1 = 5:8;
opi2 = 17:20;
figure(2)
% clf()
aa=plot(acC1(end,:), abs(sum(2*acC1(opi1,:)))/Famp, 'o-', 'LineWidth', 2); hold on
% plot(acC2(end,:), abs(sum(2*acC2(opi2,:)))/Famp, '.', 'LineWidth', 2); hold on
xlim([Wst Wen])