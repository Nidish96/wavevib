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
savdat = true;
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
Nh = length(h);
load(sprintf('./DATS/E_NLIMPDAT_2Knl%d_H%d.mat', knl,Nh), ...
     'acC', 'h', 'Famps', 'Klib', 'wcomps', ...
     'pcs', 'bcs', 'joints', 'excs', ...
     'knl', 'gap', 'Wst', 'Wen')

Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);
[zinds1,hinds1,rinds01,rinds1,iinds1] = HINDS(1, h);
Nt = 2^10;
Nh = length(h);
Nf = length(Famps);
joints.nl = @(Uw) HCONTACT(Uw, knl, gap, h, Nt);

%% Harmonic Convergence
Hmaxs = 1:101;
Nh_max = max(Hmaxs)+1;
Nhc_max = 2*max(Hmaxs)+1;
aout = zeros(Nh_max, length(Hmaxs), Nf);
ttks = zeros(length(Hmaxs), Nf);

%% Pick peak point
if analyze
    opi1 = (13:16)';
    opi2 = (17:20)';
    hsc = [sqrt(2); ones(Nh-1,1)];

    Lsel = kron(speye(Nh), sum((1:Npts*Nwc)==opi1)-sum((1:Npts*Nwc)==opi2));
    Lsel(1, [opi1 opi2]) = [1 0 0 0 -1 0 0 0];
    for fi = 1:Nf
        acout = Lsel*acC{fi}(1:end-1,:);
        [~,mi] = max(sum(hsc.*abs(acout)));

        opt = struct('ITMAX', 20, 'Display', false);
        for hi=Hmaxs
            h_n = (0:hi)';
            Nh_n = length(h_n);
            Nhc_n = sum((h_n==0)+2*(h_n~=0));
            [zinds_n,hinds_n,rinds0_n,rinds_n,iinds_n] = HINDS(Npts*Nwc, h_n);

            joints.nl = @(Uw) HCONTACT(Uw, knl, gap, h_n, Nt);

            ari0 = zeros(Npts*Nwc*Nhc_n,1);
            ari0([rinds0_n rinds_n(1:Npts*Nwc) iinds_n(1:Npts*Nwc)]) = ...
                [acC{fi}(zinds, mi); real(acC{fi}(hinds(1:Npts*Nwc), mi)); ...
                 imag(acC{fi}(hinds(1:Npts*Nwc), mi))];

            Wexc = acC{fi}(end, mi);
            tic
            aris = NSOLVE(@(ari) WVHBRESFUN([ari; Wexc], Famps(fi), h_n, pcs, bcs, joints, Klib), ...
                          ari0, opt);
            ttks(hi, fi) = toc;
            acs = zeros(Npts*Nwc*Nh_n,1);
            acs([zinds_n hinds_n]) = [aris(rinds0_n); aris(rinds_n)+1j*aris(iinds_n)];
            Lsel_n = kron(speye(Nh_n), sum((1:Npts*Nwc)==opi1)-sum((1:Npts*Nwc)==opi2));
            Lsel_n(1, [opi1 opi2]) = [1 0 0 0 -1 0 0 0];

            aout(1:Nh_n, hi, fi) = Lsel_n*acs;

            fprintf('%d: Done %d/%d\n', fi, hi, max(Hmaxs))
        end
    end

    Wexcs = zeros(1, Nf);
    for fi=1:Nf
        acout = Lsel*acC{fi}(1:end-1,:);
        [~,mi] = max(sum(hsc.*abs(acout)));

        Wexcs(fi) = acC{fi}(end, mi);
    end

    if savdat
        save('./DATS/E4_HCONV.mat', 'Hmaxs', 'Famps', 'aout', 'Wexcs')
    end
else
    load('./DATS/E4_HCONV.mat', 'Hmaxs', 'Famps', 'aout', 'Wexcs')
end

%%
fi = 2;
thresh = 1e-3;

hrelconts = squeeze(abs(aout(:,end, :)./aout(2,end, :)));
figure(1)
set(gcf, 'Color', 'white')
clf()
plot((0:Hmaxs(end)), hrelconts, '.-', 'LineWidth', 2); hold on
grid on
legend(arrayfun(@(fa) sprintf('F=%.2f N', fa), Famps, 'UniformOutput', false))
axis tight
set(gca, 'YScale', 'log')
xlabel('Harmonic Index')
ylabel('Relative Harmonic Content')

relerms = squeeze(rms(aout-aout(:,end,:))./rms(aout(:,end,:)));
hci = find(any(relerms>thresh,2), 1,'last');

figure(2)
set(gcf, 'Color', 'white')
clf()
plot(Hmaxs, relerms, '.-', 'LineWidth', 2); hold on
set(gca, 'YScale', 'log')
plot(xlim, thresh*[1 1], 'k--')
plot(Hmaxs(hci)*[1 1], ylim, 'k--')
grid on
legend(arrayfun(@(f) sprintf('F=%.2f N', f), Famps, 'UniformOutput', false), ...
       'Location', 'best')
text(Hmaxs(hci)+1, max(abs(ylim))/2, sprintf('H=%d', Hmaxs(hci)), 'fontsize', 13)

xlabel('Harmonic Truncation Order')
ylabel('RMS Harmonic Deviation')
