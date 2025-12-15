clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')
set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a portal frame with Timoshenko members
%similar to previous example.
% Here the angle joints are replaced with a spring-damper joint
%% Setup Model
Ey = 2.068e11;
nu = 0.3;
G = Ey/(2*(1+nu)); %shear modulus
kappa = 10*(1+nu)/(12+(11*nu)); %shear correction factor
rho = 7842.22747;
wid = 0.89208;
brd = 0.021255;
Ar = 0.0189612524;  % Area
Iy = 0.001257435137946;  % 2nd moment of area
%coefficients
Cb = sqrt(Ey*Iy / (rho*Ar)); %bending stiffness
Cs = sqrt(kappa*G*Ar/(rho*Ar)); %shear stiffness
Cr = sqrt(Iy / Ar); %rotational effects
% Klib K1, K2 & K3
 Klib = [struct('K', @(w,xi) sqrt( ...
     0.5*((1/Cs)^2 + (Cr/Cb)^2)*w.^2 + ...
     sqrt((w.^2)/(Cb^2) + 0.25*((1/Cs)^2 - (Cr/Cb)^2)^2 * w.^4) ));
     struct('K', @(w,xi) sqrt( ...
     abs(0.5*((1/Cs)^2 + (Cr/Cb)^2)*w.^2 - ...
     sqrt((w.^2)/(Cb^2) + 0.25*((1/Cs)^2 - (Cr/Cb)^2)^2 * w.^4))));
     struct('K', @(w,xi) w*sqrt(rho/Ey))];
%Wave components
wcomps = [-1j 1; % First component -> exp(-ik1 x )
           -1 2; % Second component-> exp(-k2 x )
           1j 1; % Third component -> exp(ik1 x )
            1 2; % Fourth component-> exp(k2 x )
          -1j 3; % Fifth component-> exp(-ik3 x )
          1j 3]; % Sixth component-> exp(ik3 x )
% The relations between the coefficients of deflection wave components & bending slope wave components
P = @(w,xi) (Klib(1).K(w,xi).^2 * Cs^2 - w.^2) ./ (Klib(1).K(w,xi) * Cs^2);
N = @(w,xi) (Klib(2).K(w,xi).^2 * Cs^2 + w.^2) ./ (Klib(2).K(w,xi) * Cs^2);
%pieces
L = 15.24;  
H = 15.24;
pcs = [struct('coords',[0 0;0 H/2;0 H],'wcomps', wcomps);
    struct('coords',[0 H;L H],'wcomps', wcomps);   
    struct('coords',[L H;L 0],'wcomps', wcomps)]; 
% BCs (Fix_Fix)
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1 0 0; -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0; 0 0 0 0 1 1]);
    struct('i', 7, 'cofs', @(w,xi) [1 1 1 1 0 0; -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0; 0 0 0 0 1 1])];
% angle-joint
m = brd*brd*wid*rho;
J = (brd^2 + brd^2)*m/12;
cofs3 = @(w,xi)[ [0 0 0 0 1 1 0 0 0 0 0 0] - [0 0 0 0 0 0 1 1 1 1 0 0];%u1-v2=0
    [0 0 0 0 0 0 0 0 0 0 1 1] + [1 1 1 1 0 0 0 0 0 0 0 0];%u2+v1=0
    [-1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0 1j*P(w,xi) 1*N(w,xi) -1j*P(w,xi) -1*N(w,xi) 0 0];% psi1-psi2=0
    (G*Ar*kappa)*[0 0 0 0 0 0 -1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0] - (Ey*Ar*Klib(3).K(w,xi))*[0 0 0 0 -1j 1j 0 0 0 0 0 0] + (m*w*w)*[0 0 0 0 1 1 0 0 0 0 0 0]; % V2-F1=0
    (G*Ar*kappa)*[-1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0 0 0 0 0 0 0] + (Ey*Ar*Klib(3).K(w,xi))*[0 0 0 0 0 0 0 0 0 0 -1j 1j] + (m*w*w)*[0 0 0 0 0 0 0 0 0 0 1 1]; % F2+V1=0
    (Ey*Iy)*[(Klib(1).K(w,xi)*P(w,xi)) (-Klib(2).K(w,xi)*N(w,xi)) (Klib(1).K(w,xi)*P(w,xi)) (-Klib(2).K(w,xi)*N(w,xi)) 0 0 (-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) (-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) 0 0] + ...
    (0.5*brd*G*Ar*kappa)*[-1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0 -1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0] +...
    (J*w*w)*[-1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0 0 0 0 0 0 0]]; % M2-M1+(h/2 * (V1+V2)) = J(w^2)psij
       
%% Non-linear spring
h = [1; 3];
Nt = 128;

kJs = diag([1e9 1e9 1e9]);
cJs = diag([5e4 5e4 5e4]);
%gJs = diag([0 0 0]); %For linear case
gJs = diag([1e14 1e14 0]);

cofs = @(w,xi) [(G*Ar*kappa)*[-1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0 0 0 0 0 0 0];
    (Ey*Iy)*[(-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) (-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) 0 0 0 0 0 0 0 0];
    Ey*Ar*Klib(3).K(w,xi)*[0 0 0 0 -1j 1j 0 0 0 0 0 0];
    (G*Ar*kappa)*[-1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0 1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) -1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 0 0];
    (Ey*Iy)*[(-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) (-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) 0 0 (Klib(1).K(w,xi)*P(w,xi)) (-Klib(2).K(w,xi)*N(w,xi)) (Klib(1).K(w,xi)*P(w,xi)) (-Klib(2).K(w,xi)*N(w,xi)) 0 0 ];
    (Ey*Ar*Klib(3).K(w,xi))*[0 0 0 0 -1j 1j 0 0 0 0 1j -1j]];

joints = [struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJs, cJs, gJs, h, Nt), ...
    'nldcofs', @(w,xi) [-1 -1 -1 -1 0 0 0 0 0 0 1 1; 0 0 0 0 -1 -1 1 1 1 1 0 0; 1j*P(w,xi) 1*N(w,xi) -1j*P(w,xi) -1*N(w,xi) 0 0 -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0], ...
    'nlfcofs', @(w,xi) [eye(3);zeros(3)]);
    struct('type', 2, 'i', 5, 'j', 6, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJs, cJs, gJs, h, Nt), ...
    'nldcofs',@(w,xi) [0 0 0 0 1 1 -1 -1 -1 -1 0 0; 1 1 1 1 0 0 0 0 0 0 -1 -1; -1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0 1j*P(w,xi) 1*N(w,xi) -1j*P(w,xi) -1*N(w,xi) 0 0], ...
    'nlfcofs', @(w,xi) [eye(3);zeros(3)])];
 

%% Excitation
Mx = @(w,xi)inv([G*Ar*kappa*[-1j*(Klib(1).K(w,xi)-P(w,xi)) (-Klib(2).K(w,xi)+N(w,xi)) 1j*(Klib(1).K(w,xi)-P(w,xi)) (Klib(2).K(w,xi)-N(w,xi)) 0 0];
                  (Ey*Iy)*[(-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) (-Klib(1).K(w,xi)*P(w,xi)) (Klib(2).K(w,xi)*N(w,xi)) 0 0];
                  Ey*Ar*Klib(3).K(w,xi)*[0 0 0 0 -1j 1j];
                  [1 1 1 1 0 0];
                  [-1j*P(w,xi) -1*N(w,xi) 1j*P(w,xi) 1*N(w,xi) 0 0];
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

Wst = 67.52;
Wen = 67.59; 
dw = 1;

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
figure(3)
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



 

    