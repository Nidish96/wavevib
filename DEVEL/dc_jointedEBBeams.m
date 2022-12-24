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
h = [1];
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
Famp = 2e3; %[1, 10, 20]

Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);

Wst = 1055.5;
Wen = 1055.9;
dw = 0.1;

Copt = struct('Nmax', 300, 'angopt', 1e-1, 'DynDscale', 1);
fmuls = [1 10 20];
acC = cell(size(fmuls));
for fi=1:length(fmuls)
    [Amat, ~, ~, Fv] = WVAMAT([Wst;0], h, pcs, bcs, joints, Klib, 'r');
    ari0 = Amat\Fv*Famp*fmuls(fi); 
    Copt.Dscale = [abs(ari0+1e-6);Wst];
    
    ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famp*fmuls(fi), h, pcs, bcs, joints, Klib), ...
        ari0, Wst, Wen, dw, Copt);
    
    acC{fi} = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
    acC{fi}([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];
end

%% 
opi = 5:8;
figure(1)
clf()
aa = gobjects(size(fmuls));
for fi=1:length(fmuls)
    aa(fi)=plot(acC{fi}(end,:), abs(sum(2*acC{fi}(opi,:)))/Famp/fmuls(fi), '-', 'LineWidth', 2); hold on
    legend(aa(fi), sprintf('F = %.2f kN', Famp*fmuls(fi)/1e3));
end
legend(aa, 'Location', 'northwest')
xlim([Wst Wen])
xlabel('Frequency (rad/s)')
ylabel('|FRF| (m/N)')

% print('./FIGS/fresp_EBBEAMJ.png', '-dpng')

%% Evaluate Mode Shape
Xs = linspace(0, 2*L0, 1000)';
dL = [1 1 1 1];

figure(2)
clf()
for fi=1:length(fmuls)

    [~, iw] = max(abs(sum(2*acC{fi}(opi,:))));
    Vs = WVEVALWCOFS(acC{fi}(1:Npts*Nwc,iw)*exp(-1j*angle(sum(acC{fi}(opi,iw)))), acC{fi}(end,iw), dL, Xs, pcs, Klib);
    
    subplot(2,2,fi)
    plot(Xs, real(Vs)/Famp/fmuls(fi), 'LineWidth', 2); hold on
    plot(Xs, imag(Vs)/Famp/fmuls(fi), 'LineWidth', 2)
    xlabel('X Coordinate (m)')
    ylabel('v/F (m/N)')
    for i=1:length(pcs)
        k1 = pcs(i).irange(1);
        k2 = pcs(i).irange(2);
    
        AC = kron(eye(k2-k1+1), dL)*acC{fi}((k1-1)*Nwc+1:k2*Nwc,iw)*exp(-1j*angle(sum(acC{fi}(opi,iw))));
        plot(pcs(i).coords(:,1), real(AC)/Famp/fmuls(fi), 'ko', 'MarkerFaceColor', 'w');
        plot(pcs(i).coords(:,1), imag(AC)/Famp/fmuls(fi), 'ko', 'MarkerFaceColor', 'w');
    end
    title(sprintf('%.2f kN', Famp*fmuls(fi)/1e3))
end

% print('./FIGS/dshape_EBBEAMJ.png', '-dpng')