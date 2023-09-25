% Conduct harmonic convergence of wave-based model
clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')

addpath('../ROUTINES/export_fig')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a nonlinear jointed Euler-Bernoulli Beam example
%similar to the jointed bar example. Parameters taken from 
% Krishna and Chandramouli, 2012

savfig = false;
animfig = false;
savdat = false;
analyze = false;
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
knl = 1210/2;      % 220,      880,    (242, 484, 880, 1210)
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
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j], ...
              'cofs0', @(xi) [eye(2) zeros(2)]);  % fixed end
       struct('i', 5, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j], ...
              'cofs0', @(xi) [eye(2) zeros(2)]);     % fixed end
       struct('i', 7, 'cofs', @(w,xi) [1 1 1 1; 1 -1 1j -1j], ...
              'cofs0', @(xi) [eye(2) zeros(2)]);     % fixed end
       struct('i', 3, 'cofs', @(w,xi) [1 1 -1 -1], ...
              'cofs0', @(xi) [0 0 1 0]);        % Moment-free
       struct('i', 4, 'cofs', @(w,xi) [1 1 -1 -1], ...
              'cofs0', @(xi) [0 0 1 0]);        % Moment-free
       struct('i', 6, 'cofs', @(w,xi) [1 1 -1 -1], ...
              'cofs0', @(xi) [0 0 1 0])];       % Moment-free

%% Setup Excitation
Mx = @(w,xi) inv([Ey*Iy*Klib.K(w,xi)^3*[-1 1 1j -1j];
                  Ey*Iy*Klib.K(w,xi)^2*[-1 -1 1 1];
                  [1 1 1 1];
                  Klib.K(w,xi)*[1 -1 1j -1j]]);
excs = struct('i', 2, 'nh', 1, ...
              'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0], ...
              'rcofs0', @(xi) [0;0;0;1/(6*Ey*Iy)]);
%'nh' sets the harmonic at which to apply the excitation

%% Setup Joints
cofs = @(w,xi) [-[zeros(1,4) 1 -1 -1j 1j zeros(1,4)];
    -[zeros(1,8) 1 -1 -1j 1j];
    kron([1 -1 -1], [1 -1 -1j 1j])];
joints = struct('type', 3, 'is', [3 4 6], ...
                'cofs', cofs, 'nl', [], ...  % Setup nl after settin 'h'
    'nldcofs', @(w,xi) kron([1 -1 0;1 0 -1], [1 1 1 1]), ...
    'nlfcofs', @(w,xi) [eye(2);0 0]/(Ey*Iy*Klib.K(w,xi)^3));

%% Pre-Processing
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%% Conduct Nonlinear Forced Response Analysis
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

h = (0:16)';
Nt = 2^10;
Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);
[zinds1,hinds1,rinds01,rinds1,iinds1] = HINDS(1, h);

load(sprintf('./DATS/E_NLIMPDAT_2Knl%d_H%d.mat', knl,Nh), ...
     'acC', 'h', 'Famps', 'Klib', 'wcomps', ...
     'pcs', 'bcs', 'joints', 'excs', ...
     'knl', 'gap', 'Wst', 'Wen')

Nf = length(Famps);
joints.nl = @(Uw) HCONTACT(Uw, knl, gap, h, Nt);

%% Load FE Results
Nn1 = 10;
fe = load(sprintf('./DATS/E1_FERES_2Knl%d_N%d_H%d.mat', 2*knl, Nn1, max(h)), ...
          'UCs', 'Fl', 'Famps', 'Rb', 'h', 'Nt', 'Wst', 'Wen', 'Nn1', 'Lb');
fe.Rrel = fe.Lb((fe.Nn1-1)*2+1,:)-fe.Lb((fe.Nn1+1-1)*2+1,:);

%% Plot Forced Response
opi1 = (13:16)';  % x1
opi2 = (17:20)';  % x2
opi3 = (25:28)';  % x3
hi = find(h==1);  % Choose which harmonic to plot

% fe.Rout = fe.Rb;
% wbm.Rout = 2*sum((1:Npts*Nwc)==opi1);  % Abs response at x1

fe.Rout = fe.Rrel;
wbm.Rout = 2*(sum((1:Npts*Nwc)==opi1)-sum((1:Npts*Nwc)==opi2));  % relative coords

colos = DISTINGUISHABLE_COLORS(length(Famps))*0.75;

figure(1)
set(gcf, 'Color', 'white')
clf()
aa = gobjects(size(fe.Famps));
for fi=1:length(fe.Famps)
    if fi<=length(fe.Famps)
        Uh = kron(eye(Nhc), fe.Rout)*fe.UCs{fi}(1:end-1,:);
    end
    aout = wbm.Rout*acC{fi}((hi-1)*Npts*Nwc+(1:Npts*Nwc),:);
    
    subplot(2,1,1)
    aa(fi)=plot(acC{fi}(end,:)/2/pi, abs(aout)*1e3, ...
                '-', 'LineWidth', 2, 'Color', colos(fi,:)); hold on
    legend(aa(fi), sprintf('F = %.1f N', Famps(fi)));
    if fi<=length(fe.Famps)
        bb = plot(fe.UCs{fi}(end,:)/2/pi, abs([1 1j]*Uh(2:3,:)*1e3), ...
                  '--', 'LineWidth', 2, 'Color', colos(fi,:))
        legend(bb, 'FE')
    end
    grid on
    ylabel(sprintf('H%d Response (mm)', h(hi)))
    subplot(2,1,2)
    plot(acC{fi}(end,:)/2/pi, rad2deg(angle(aout)), ...
         '-', 'LineWidth', 1.5); hold on
    if fi<=length(fe.Famps)
        plot(fe.UCs{fi}(end,:)/2/pi, rad2deg(angle([1 1j]*Uh(2:3,:))), ...
             '--', 'LineWidth', 2, 'Color', colos(fi,:))
    end
    grid on
    ylabel('Phase (degs)')
end
subplot(2,1,1)
if h(hi)==1
    plot(xlim, gap*[1 1]*1e3, 'k-.')
    text(10, gap*1.25*1e3, 'Gap limit')
end
xlim(sort([Wst Wen]/2/pi))
legend([aa bb], 'Location', 'northwest')
subplot(2,1,2)
xlim(sort([Wst Wen]/2/pi))
set(gca, 'YTick', -180:45:180)
ylim([0 1]*180)
xlabel('Frequency (Hz)')

%% Plot the response at one point (From WBM Solution)
t = linspace(0, 2*pi, Nt+1)'; t=t(1:end-1);
fi = 2;
opis = {(13:16), (17:20), (25:28)};  % Points x1, x2, x3 in paper

Uouts = cell(size(opis));
uouts = cell(size(opis));
for i=1:length(opis)
    Rout = speye(Npts*Nwc);
    Rout = kron(eye(Nh), 2*sum(Rout(opis{i},:)));
    if h(1)==0
        Rout(1, opis{i}) = [1 0 0 0];
    end
    
    Uouts{i} = Rout*acC{fi}(1:end-1,:);
    Uouth = zeros(Nhc, size(Uouts{i},2));
    Uouth([rinds01 rinds1 iinds1], :) = ...
        [real(Uouts{i}(zinds1,:));
         real(Uouts{i}(hinds1,:)); imag(Uouts{i}(hinds1,:))];
    uouts{i} = AFT(Uouth, h,Nt, 'f2t');
    if i>1
        uouts{i} = uouts{i}+(-1)^i*gap;
    end
end
[~, mi] = max(max(abs(uouts{1})));
% [~, mi] = min(abs(acC{fi}(end,:)/2/pi-9));

colos = DISTINGUISHABLE_COLORS(length(opis))*0.75;
aa = gobjects(size(opis));
figure(2)
set(gcf, 'Color', 'white')
clf()
subplot(2,1,1)
for i=1:length(opis)
    aa(i) = plot(acC{fi}(end,:)/2/pi, max(abs(uouts{i}))*1e3, '-', ...
                 'LineWidth', 1.5, 'Color', colos(i,:)); hold on
    plot(acC{fi}(end,mi)/2/pi, max(abs(uouts{i}(:,mi)))*1e3, 'ko', ...
         'MarkerFaceColor', 'k');

    legend(aa(i), ['Point x' num2str(i)])
end
legend(aa, 'Location', 'northwest')
axis tight
grid on
xlabel('Frequency (Hz)')
ylabel('Peak Amplitude (mm)')
subplot(2,1,2)
for i=1:length(opis)
    plot(t, uouts{i}(:, mi)*1e3, 'LineWidth', 2, ...
         'Color', colos(i,:)); hold on
end
set(gca, 'XTick', 0:pi/2:2*pi, 'XTicklabel', ...
         {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
axis tight
xlim([0 2*pi])
plot(xlim, gap*[1 1]*1e3, 'k--'); hold on
plot(xlim, -gap*[1 1]*1e3, 'k--'); hold on
grid on
xlabel('Scaled Time')
ylabel('Response')
