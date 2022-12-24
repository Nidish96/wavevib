clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This shows the non-linear forced response computation for 
%the 1D bar example considered in Balaji, Brake, Leamy (2022a,2022b)

%% Setup model
Ey = 2.62e11;
rho = 1280;
ell = 1.0;
Ar = 31.75e-3*54e-3;
Klib = struct('K', @(w,xi) w/sqrt(Ey/rho));
wcomps = [1j 1;  % First component -> exp( j k x )
         -1j 1]; % Sec component   -> exp(-j k x )

% Setup "wave-based pieces"
pcs = [struct('coords', [0;0.28;ell/3], 'wcomps', wcomps);
    struct('coords', [ell/3;ell], 'wcomps', wcomps)];

% Setup Boundary Conditions (pts indexed as given above in sequence)
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1]);
    struct('i', 5, 'cofs', @(w,xi) [1 1])];

%% Setup the Nonlinear Joint
h = sort([1 3]);  % List of harmonics to consider for simulation. 
Nt = 128;  % Number of points for Alternating Frequency-Time calculations.

kJ = 1e9;
cJ = 320;
gJ = 1e8;
cofs = @(w,xi) [1j*Klib.K(w,xi)*[1 -1 0 0];1 -1 -1 1];
joints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJ, cJ, gJ, h, Nt), ...
    'nldcofs', @(w,xi) [1 1 -1 -1], ... % The relative displacement.
    'nlfcofs', @(w,xi) 1/(Ey*Ar)*[1; 0]);
%NOTE: The 'nldcofs' specifies the coefficient matrix such that the
% displacement(s) used as input for the nonlinearity is given as,
%   u = NLDCOFS * [a;b]         where "a" and "b" are the wave coefficients
%   of points 'i' and 'j' respectively.
%
%      The 'nlfcofs' specifies the coefficient matrix that the results from
% the nonlinear force function are multiplied with before being added to
% the equations of motion. The nonlinear forces are added in the LHS of the
% equations. 

%% Setup Excitation
excs = struct('i', 2, ...
    'rcofs', @(w,xi) (1/2/(2j*Klib.K(w,xi)*Ey*Ar))*[-1;1]);

%% Preprocess Everything
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);
Nwc = size(wcomps,1);  % Number of wave components

%% Conduct Nonlinear Forced Response Analysis using Continuation
Npts = pcs(end).irange(end);
Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
% Compute indices to convert between real-imaginary and complex
% representations:
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);

Wst = 78e3;
Wen = 82e3;
dw  = 0.025;

Copt = struct('Nmax', 100, 'angopt', 2e-1, 'DynDscale', 1);

Famps = 150e5*[0.5 1 2];
acC = cell(size(Famps));
for fi=1:length(Famps)
    % Continuation is done in real-imaginary representation
    ari0 = zeros(Npts*Nwc*Nhc, 1);
    ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famps(fi), h, pcs, bcs, joints, Klib), ...
        ari0, Wst, Wen, dw, Copt);
    
    % Convert to complex representation
    acC{fi} = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
    acC{fi}([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];
end

%%
opi =  9:10;
figure(2)
clf()
for fi=1:length(Famps)
    subplot(2,1,1)
    aa(fi)=plot(acC{fi}(end,:)/1e3, abs(2*sum(acC{fi}(opi,:))), 'LineWidth', 2); hold on
    legend(aa(fi), sprintf('F = %.1f MN', Famps(fi)*1e-6))
    ylabel('Response (m)')
    subplot(2,1,2)
    plot(acC{fi}(end,:)/1e3, rad2deg(angle(2*sum(acC{fi}(opi,:)))), 'LineWidth', 2); hold on
    ylabel('Phase (degs)')
end
for i=1:2
    subplot(2,1,i)
    xlim([79 82]); grid on
end
set(gca, 'YTick', -180:90:180)
xlabel('Frequency (k rad/s)')
subplot(2,1,1)
legend(aa, 'Location', 'northwest')
