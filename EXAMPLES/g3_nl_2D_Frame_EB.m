clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a portal frame with EB members example 
% similar to previous examples.
% Here the angle joints are replaced with nonlinear spring-damper joints.

%% Setup Model
Ey = 2.068e11;
rho = 7842.22747;
wid = 0.89208;
brd = 0.021255;
Ar = 0.0189612524;  % Area
Iy = 0.001257435137946;  % 2nd moment of area
Klib = [struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25));
        struct('K', @(w,xi) w*sqrt(rho/Ey))];
%Wave components
wcomps = [1 1;  % -> exp(k1 x )
         -1 1;  % -> exp(-k1 x )
         1j 1;  % -> exp(ik1 x )
        -1j 1;  % -> exp(-ik1 x )
         1j 2;  % -> exp(ik2 x )
        -1j 2]; % -> exp(-ik2 x )
%pieces
L = 15.24;  
H = 15.24; 
pcs = [struct('coords',[0 0;0 H/2;0 H],'wcomps', wcomps);
       struct('coords',[0 H;L H],'wcomps', wcomps);   
       struct('coords',[L H;L 0],'wcomps', wcomps)]; 
% BCs (Fix_Fix)
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1 0 0; 1 -1 1j -1j 0 0; 0 0 0 0 1 1]);
       struct('i', 7, 'cofs', @(w,xi) [1 1 1 1 0 0; 1 -1 1j -1j 0 0; 0 0 0 0 1 1])];

%% Non-linear spring (instead of rigid angle joints)
h = [1; 3];
Nt = 128;
kJs = diag([1e9 1e9 1e9]);
cJs = diag([5e4 5e4 5e4]);
gJs = diag([5e15 5e15 0]);

cofs = @(w,xi) [(-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0 0 0 0 0 0 0];
    (Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0 0 0 0 0 0 0];
    (Ey*Ar*Klib(2).K(w,xi))*[0 0 0 0 1j -1j 0 0 0 0 0 0];
    (-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0 -1 1 1j -1j  0 0];
    (Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0 -1 -1 1 1 0 0];
    (Ey*Ar*Klib(2).K(w,xi))*[0 0 0 0 1j -1j 0 0 0 0 -1j 1j]];

joints = [struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJs, cJs, gJs, h, Nt), ...
    'nldcofs', @(w,xi) [-1 -1 -1 -1 0 0 0 0 0 0 1 1; ...
    0 0 0 0 -1 -1 1 1 1 1 0 0; Klib(1).K(w,xi)*[-1 1 -1j 1j 0 0 1 -1 1j -1j 0 0]], ...
    'nlfcofs', @(w,xi) [eye(3);zeros(3)]);
    struct('type', 2, 'i', 5, 'j', 6, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJs, cJs, gJs, h, Nt), ...
    'nldcofs', @(w,xi) [0 0 0 0 1 1 -1 -1 -1 -1 0 0; ...
    1 1 1 1 0 0 0 0 0 0 -1 -1; Klib(1).K(w,xi)*[1 -1 1j -1j 0 0 -1 1 -1j 1j 0 0]], ...
    'nlfcofs', @(w,xi) [eye(3);zeros(3)])];

%% Excitation
Mx = @(w,xi)inv([-(Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0];
    (Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0];
    (Ey*Ar*Klib(2).K(w,xi))*[0 0 0 0 1j -1j];
    [1 1 1 1 0 0];
    Klib(1).K(w,xi)*[1 -1 1j -1j 0 0];
    [0 0 0 0 1 1]]);

excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0;0;0], ...
    'rcofs0', [1/2;0;0;0;0;0]);

%% Pre-Processing
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%% Conduct Nonlinear Forced Response Analysis using Continuation
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);

Wst = 65.5; 
Wen = 65.9; 
dw = 0.5; 

Copt = struct('Nmax', 700, 'angopt', 1e-2, 'DynDscale', 1);
Famps = 1e2*[1 5 10];
acC = cell(size(Famps));
for fi=1:length(Famps)
    ari0 = zeros(Npts*Nwc*Nhc, 1);
    Copt.Dscale = [abs(ari0+1e-2);Wst];  % Setup Scaling for continuation (sometimes helps)

    % Conduct continuation
    ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famps(fi), h, pcs, bcs, joints, Klib), ...
        ari0, Wst, Wen, dw, Copt);

    % Convert to complex representation
    acC{fi} = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
    acC{fi}([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];
end

%% Plot Results
opi = 7:12;
figure()
clf()
aa = gobjects(size(Famps));
for fi=1:length(Famps)
    subplot(2,1,1)
    aa(fi)=plot(acC{fi}(end,:), abs(sum(2*acC{fi}(opi,:))), '-', 'LineWidth', 2); hold on
    legend(aa(fi), sprintf('F = %.0f N', Famps(fi)));
    grid on
    ylabel('Response (m)')
    subplot(2,1,2)
    plot(acC{fi}(end,:), rad2deg(angle(sum(2*acC{fi}(opi,:)))), '-', 'LineWidth', 2); hold on
    grid on
    ylabel('Phase (degs)')
end
subplot(2,1,1)
xlim([Wst Wen])
legend(aa, 'Location', 'northwest')
subplot(2,1,2)
xlim([Wst Wen])
set(gca, 'YTick', -180:90:180)
xlabel('Frequency (rad/s)')



 

    