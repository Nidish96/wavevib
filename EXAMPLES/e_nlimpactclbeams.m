clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a nonlinear jointed Euler-Bernoulli Beam example
%similar to the jointed bar example. Parameters taken from 
% Krishna and Chandramouli, 2012

savfig = true;
%% Setup Model
Ey = 2.1e11;
rho = 7680;
thk = 3e-3;   % Thickness 1
wid  = 40e-3;  % Width
Ar = thk*wid;  % Area
Iy = thk^3*wid/12;  % 2nd moment of area
L1  = 560e-3;  % Primary beam length
xF1 = L1/6;    % Excitation Location (on primary beam)
L2 = 445e-3;  % Secondary beam length (each)
% Klib = struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25));  % Undamped EB-Beam Case
al = 0.80;    % 0.80, 1.80
bt = 1.1e-4;  % 1.1e-4, 2.475e-4
% Joint parameters
knl = 220/2;      % 220,      880,    (242, 484, 880, 1210)
gap = 2.5e-3;   % 2.5e-3,   0.35,   0.35

Klib = struct('K', @(w,xi) ((rho*Ar*(w.^2+1j*w*al))./(Ey*Iy*(1-1j*w*bt))).^(0.25) );
wcomps = [1 1;  % First component -> exp(  k x )
         -1 1;  % Second component-> exp( -k x )
         1j 1;  % Third component -> exp( ik x )
        -1j 1]; % Fourth component-> exp(-ik x )

% Setup "wave-based pieces"
pcs = [struct('coords', [0 0;xF1 0;L1 0], 'wcomps', wcomps);
    struct('coords', [L1 gap;L1+L2 gap], 'wcomps', wcomps);
    struct('coords', [L1 -gap;L1+L2 -gap], 'wcomps', wcomps)];

% Setup Boundary Conditions. Fix-Fix used here.
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j]);  % fixed end
    struct('i', 5, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j]);     % fixed end
    struct('i', 7, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j]);     % fixed end
    struct('i', 3, 'cofs', @(w,xi) [1 1 -1 -1]);        % Moment-free
    struct('i', 4, 'cofs', @(w,xi) [1 1 -1 -1]);        % Moment-free
    struct('i', 6, 'cofs', @(w,xi) [1 1 -1 -1])];       % Moment-free

%% Setup Excitation
Mx = @(w,xi) inv([Ey*Iy*Klib.K(w,xi)^3*[-1 1 1j -1j];
                  Ey*Iy*Klib.K(w,xi)^2*[-1 -1 1 1];
                  [1 1 1 1];
                  Klib.K(w,xi)*[1 -1 1j -1j]]);

excs = struct('i', 2, 'nh', 1, ...
    'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0]);
%'nh' sets the harmonic at which to apply the excitation

%% Setup AFT parameters
h = (0)';
Nt = 2^10;

%% Setup Joints
cofs = @(w,xi) [-[zeros(1,4) 1 -1 -1j 1j zeros(1,4)];
    -[zeros(1,8) 1 -1 -1j 1j];
    kron([1 -1 -1], [1 -1 -1j 1j])];
joints = struct('type', 3, 'is', [3 4 6], ...
    'cofs', cofs, 'nl', @(Uw) HCONTACT(Uw, knl, gap, h, Nt), ...
    'nldcofs', @(w,xi) kron([1 -1 0;1 0 -1], [1 1 1 1]), ...
    'nlfcofs', @(w,xi) [eye(2);0 0]/(Ey*Iy*Klib.K(w,xi)^3));

%% Pre-Processing
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%% Conduct Nonlinear Forced Response Analysis using Continuation
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);

Wst = 2*pi*4;
Wen = 2*pi*12;
dw = 0.75;

Copt = struct('Nmax', 600, 'angopt', 5e-2, 'DynDscale', 1);
Famps = [0.4 1.0 2.5];  % 0.336, 0.336
acC = cell(size(Famps));
for fi=1:length(Famps)
    [Amat, ~, ~, Fv, ~, ~, JEV] = WVAMAT([Wst;0], h, pcs, bcs, joints, Klib, 'r');
    
    ari0 = Amat\(Fv*Famps(fi));

    Copt.Dscale = [1e-6*ones(size(ari0));Wst];
    ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famps(fi), h, pcs, bcs, joints, Klib), ...
        ari0, Wst, Wen, dw, Copt);

    % Convert to complex representation
    acC{fi} = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
    acC{fi}([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];
end

%% Plot Results
opi = (13:16);
hi = find(h==1);  % Choose which harmonic to plot
figure(1)
set(gcf, 'Color', 'white')
clf()
aa = gobjects(size(Famps));
for fi=1:length(Famps)
    subplot(2,1,1)
    aa(fi)=plot(acC{fi}(end,:)/2/pi, abs(sum(2*acC{fi}((hi-1)*Npts*4+opi,:)))*1e3, '-', 'LineWidth', 1.5); hold on
    legend(aa(fi), sprintf('F = %.1f N', Famps(fi)));
    grid on
    ylabel(sprintf('H%d Response (mm)', h(hi)))
    subplot(2,1,2)
    plot(acC{fi}(end,:)/2/pi, rad2deg(angle(sum(2*acC{fi}((hi-1)*Npts*4+opi,:)))), '-', 'LineWidth', 1.5); hold on
    grid on
    ylabel('Phase (degs)')
end
subplot(2,1,1)
if h(hi)==1
    plot(xlim, gap*[1 1]*1e3, 'k-.')
    text(10, gap*1.25*1e3, 'Gap limit')
end
xlim(sort([Wst Wen]/2/pi))
legend(aa, 'Location', 'northwest')
subplot(2,1,2)
xlim(sort([Wst Wen]/2/pi))
set(gca, 'YTick', -180:45:180)
ylim([0 1]*180)
xlabel('Frequency (Hz)')

if savfig
    print('./FIGS/E_out.png', '-dpng', '-r300')
end

%% Plot deflection shape
fip = 3;
wip = 180;

ac = acC{fip}(1:Npts*4, wip);
w = acC{fip}(end, wip);

figure(2)
clf()
cbk(fip,1,Famps,acC,Npts,opi,sort([Wst Wen]/2/pi),pcs,Klib,h);
sl = uicontrol('Style', 'slider', 'Position', [85 425 425 15], ...
    'String', 'Point', 'Value', 1, 'Min', 1, 'max', size(acC{fip},2), ...
    'SliderStep', 1/size(acC{fip},2)*[1 1], ...
    'Callback', @(es,ed) cbk(fip,fix(es.Value),Famps,acC,Npts,opi,sort([Wst Wen]/2/pi),pcs,Klib,h));

function [] = cbk(fip,wip,Famps,acC,Npts,opi,Wrng,pcs,Klib,h)
    ac = reshape(acC{fip}(1:end-1, wip), Npts*4, []);
    w = acC{fip}(end, wip);
    
    Nx = 10;
    [Us, Xs] = WVEVALWCOFS(ac, w, h, ones(1,4), Nx, pcs, Klib);
    
    uall = cell2mat(Us);
    uall = uall(:);
%     uall = uall./uall;
    [~,mi] = max(abs(uall));
    Us = cellfun(@(u) u*exp(-1j*angle(uall(mi))), Us, 'UniformOutput', false);
    
    hi = 1;
    colos = eye(3)*0.8;
    subplot(2,1,1)
    cla()
    aa = gobjects(size(Famps));
    for fi=1:length(Famps)
        aa(fi)=plot(acC{fi}(end,:)/2/pi, abs(sum(2*acC{fi}((hi-1)*Npts*4+opi,:)))*1e3, '-', 'LineWidth', 1.5); hold on
        legend(aa(fi), sprintf('F = %.1f N', Famps(fi)));
    end
    plot(acC{fip}(end,wip)/2/pi, abs(sum(2*acC{fip}((hi-1)*Npts*4+opi,wip)))*1e3, 'ko', 'MarkerFaceColor', 'k'); hold on
    xlim(Wrng)
    grid on
    ylabel(sprintf('H%d Response (mm)', h(hi)))
    xlabel('Frequency (Hz)')
    legend(aa, 'Location', 'northwest')
    subplot(2,1,2)
    cla()
    for n=1:length(pcs)
        plot(Xs{n}(:,1), Xs{n}(:,2), 'k-'); hold on
    
        plot(Xs{n}(:,1), Xs{n}(:,2)+real(Us{n}(:,hi)), '.-', 'Color', colos(n,:))
        plot(Xs{n}(:,1), Xs{n}(:,2)+imag(Us{n}(:,hi)), '.--', 'Color', colos(n,:))
    
        achi = ac(:, hi);
        Up = sum(achi(((pcs(n).irange(1):pcs(n).irange(2))'-1)*4+(1:4)), 2)*exp(-1j*angle(uall(mi)));
        plot(pcs(n).coords(:,1), pcs(n).coords(:,2)+imag(Up), 'ro', 'MarkerFaceColor', 'w')
        plot(pcs(n).coords(:,1), pcs(n).coords(:,2)+real(Up), 'ko', 'MarkerFaceColor', 'r')
    end
    grid on
    xlim([0 pcs(3).coords(end,1)])
    xlabel('X Coordinate (m)')
    ylabel('Y Coordinate (m)')
end
