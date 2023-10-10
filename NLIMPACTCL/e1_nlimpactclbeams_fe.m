clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')
addpath('../ROUTINES/FEM/BEAMS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%DESCRIPTION: This is a nonlinear jointed Euler-Bernoulli Beam example
%similar to the jointed bar example. Parameters taken from 
% Krishna and Chandramouli, 2012

animfig = false;
savfig = false;
savdat = true;
analyze = true;

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
knl = 1210/2;    % 220,      880,    (242, 484, 880, 1210)
gap = 2.5e-3;   % 2.5e-3,   0.35,   0.35

Nn1 = 20;
Nn2 = fix(.8*Nn1);
Xn1 = [linspace(0, L1, Nn1)' zeros(Nn1,1)];
Xn2 = [linspace(L1, L1+L2, Nn2)' gap*ones(Nn2,1)];
Xn3 = [linspace(L1, L1+L2, Nn2)' -gap*ones(Nn2,1)];

[~,mi]=min(abs(Xn1(:,1)-xF1));
Xn1(mi,1) = xF1;
Xn1(:,1) = sort(Xn1(:,1));
iF1 = find(Xn1(:,1)==xF1);

%% Construct Matrices
M1 = zeros(Nn1*2);
K1 = zeros(Nn1*2);
for e=1:Nn1-1
    [Me, Ke] = EBBEAM_MATS(rho, Ey, Ar, Iy, Xn1(e+1,1)-Xn1(e,1));

    is = (e-1)*2+1;
    ie = (e+1)*2;
    M1(is:ie, is:ie) = M1(is:ie, is:ie) + Me([2 3 5 6], [2 3 5 6]);
    K1(is:ie, is:ie) = K1(is:ie, is:ie) + Ke([2 3 5 6], [2 3 5 6]);
end
M1 = sparse(M1);
K1 = sparse(K1);

M2 = zeros(Nn2*2);
K2 = zeros(Nn2*2);
for e=1:Nn2-1
    [Me, Ke] = EBBEAM_MATS(rho, Ey, Ar, Iy, Xn2(e+1,1)-Xn2(e,1));

    is = (e-1)*2+1;
    ie = (e+1)*2;
    M2(is:ie, is:ie) = M2(is:ie, is:ie) + Me([2 3 5 6], [2 3 5 6]);
    K2(is:ie, is:ie) = K2(is:ie, is:ie) + Ke([2 3 5 6], [2 3 5 6]);
end
M2 = sparse(M2);
K2 = sparse(K2);

%% Construct full System
Mf = blkdiag(M1, M2, M2);
Kf = blkdiag(K1, K2, K2);

% BCs
Lb = eye(Nn1+2*Nn2);
Lb(:, [1 Nn1+Nn2 Nn1+2*Nn2]) = [];
Lb = sparse(kron(Lb, eye(2)));

Mb = Lb'*Mf*Lb;
Kb = Lb'*Kf*Lb;
Cb = al*Mb + bt*Kb;  % Check this!
Fv = Lb((iF1-1)*2+1,:)';
Rb = Lb((Nn1-1)*2+1,:);  % Response pt -> last node of primary beam

%% Setup Nonlinearity
nl1i = ([Nn1+1 Nn1]-1)*2+1;
nl2i = ([Nn1 Nn1+Nn2+1]-1)*2+1;

Lnls = [diff(Lb(nl1i,:));
    diff(Lb(nl2i,:))];

nlfunc = @(t, u, ud) deal(max(knl*(u-gap),0), knl*(u>gap), zeros(size(u)));

MDL = MDOFGEN(Mb,Kb,Cb,Lb);
MDL = MDL.SETNLFUN(1+3, Lnls, nlfunc);

%%
h = (0:3)';  % 16 is "converged"
Nt = 1024;
Nhc = sum((h==0)+2*(h~=0));
[~,~,zinds,rinds,iinds] = HINDS(1,h);

Fl = zeros(Nhc,1);
Fl(rinds(1)) = 1;
Fl = kron(Fl, Fv);

Famps = [0.4 1.0 2.5];  % 0.336, 0.336

Wst = 2*pi*4;
Wen = 2*pi*12;
dw  = 0.75;

Copt = struct('Nmax', 5000, 'angopt', 1e-2, 'DynDscale', 1, ...
    'Dscale', [1e-5*ones(MDL.Ndofs*Nhc,1);Wst]);

if analyze
    UCs = cell(size(Famps));
    for fi=1:length(Famps)
        U0 = HARMONICSTIFFNESS(MDL.M,MDL.C,MDL.K,Wst,h)\(Fl*Famps(fi));

        UCs{fi} = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl*Famps(fi), h, Nt), ...
                           U0, Wst, Wen, dw, Copt);
    end
    if savdat
        save(sprintf('./DATS/E1_FERES_2Knl%d_N%d_H%d.mat', 2*knl, Nn1, max(h)), ...
             'UCs', 'Fl', 'Famps', 'Rb', 'h', 'Nt', 'Wst', 'Wen', 'Nn1', 'Lb');
    end
else
    load(sprintf('./DATS/E1_FERES_2Knl%d_N%d_H%d.mat', 2*knl, Nn1, max(h)), ...
         'UCs', 'Fl', 'Famps', 'Rb', 'h', 'Nt', 'Wst', 'Wen', 'Nn1', 'Lb')
end

%% Plotting
figure(1)
clf()
set(gcf, 'Color', 'white')
aa = gobjects(size(Famps));
for fi=1:length(Famps)
    Uh = kron(eye(Nhc), Rb)*UCs{fi}(1:end-1,:);

    subplot(2,1,1)
    aa(fi) = plot(UCs{fi}(end,:)/2/pi, abs([1 1j]*Uh(2:3,:)*1e3), '-', 'LineWidth', 2); hold on
    legend(aa(fi), sprintf('F = %.1f N', Famps(fi)))

    subplot(2,1,2)
    plot(UCs{fi}(end,:)/2/pi, rad2deg(angle([1 1j]*Uh(2:3,:))), '-', 'LineWidth', 2); hold on
end
subplot(2,1,1)
plot(xlim, gap*1e3*[1 1], 'k-.')
text(10, gap*1.25*1e3, 'Gap limit')
ylabel('H1 Response (mm)')
xlim([Wst Wen]/2/pi)
legend(aa, 'Location', 'northwest')
grid on

subplot(2,1,2)
xlabel('Frequency (Hz)')
ylabel('Phase (degs)')
xlim([Wst Wen]/2/pi)
grid on
set(gca, 'YTick', -180:45:180)
ylim([0 1]*180)

if savfig
    print('./FIGS/E_feout.png', '-dpng', '-r300')

    savefig('./FIGS/E_fefig.fig')
end

%% Higher Harmonics Figure
Nh = max(h);
% Rout = Rb;
Rout = Lnls(1, :);

figure(10)
pos = get(gcf, 'Position');
set(gcf, 'Color', 'white')
clf()
aa = gobjects(size(Famps));
subplot(6, 3, 1)
for fi=1:length(Famps)
    Uh = kron(eye(Nhc), Rout)*UCs{fi}(1:end-1,:)*1e3;

    aa(fi) = plot(UCs{fi}(end,:)/2/pi, sum(abs(Uh))/sqrt(2), ...
                  '-', 'LineWidth', 2); hold on
    legend(aa(fi), sprintf('F = %.1f N', Famps(fi)))
end
plot(xlim, [1 1]*gap*1e3/sqrt(2), 'k--');
legend(aa, 'Location', 'northwest');
xlim([Wst Wen]/2/pi)
grid on
ylabel('RMS')
title('Interference Response (mm)')

hord = reshape(reshape(0:length(h), 6,3)', [],1);
[~, hord] = sort(hord);
for hi=1:length(h)
    subplot(6,3, hord(1+hi))
    for fi=1:length(Famps)
        Uh = kron(eye(Nhc), Rout)*UCs{fi}(1:end-1,:)*1e3;

        inds = 1+(h(hi)-1)*2+(1:2);
        inds = inds(inds>0);
        wv = [1 1j];
        if h(hi)==0
            wv = 1;
        end
        plot(UCs{fi}(end,:)/2/pi, abs(wv*Uh(inds, :)), ...
             '-', 'LineWidth', 2); hold on
    end
    if hi==1
        plot(xlim, [1 1]*gap*1e3, 'k--');
    end

    xlim([Wst Wen]/2/pi)
    ylabel(sprintf('H%d', h(hi)))
    if hi~=Nh
        set(gca, 'XTickLabel', [])
    end
    grid on
end
xlabel('Frequency (Hz)')
if savfig
    print('./FIGS/E_feout_harms.png', '-dpng', '-r300')

    savefig('./FIGS/E_fefig_harms.fig')
end

%% Time domain @ peak
t = linspace(0, 2*pi, Nt+1)'; t=t(1:end-1);
fi = 2;

ns = ([Nn1 Nn1+1 Nn1+Nn2+1]-1)*2+1;
uouts = cell(size(ns));
for i=1:length(ns)
    uouts{i} = AFT(kron(eye(Nhc), Lb(ns(i), :))*UCs{fi}(1:end-1,:), h,Nt, 'f2t');
    if i>1
        uouts{i} = uouts{i} + (-1)^(i)*gap;
    end
end
[~, mi] = max(max(abs(uouts{1})));
% [~, mi] = min(abs(UCs{fi}(end,:)/2/pi-10));
mi = mi;

figure(100)
clf()
subplot(2,1,1)
plot(UCs{fi}(end,:)/2/pi, max(abs(uouts{1}))*1e3, ...
     'LineWidth', 2); hold on
plot(UCs{fi}(end,mi)/2/pi, max(abs(uouts{1}(:,mi)))*1e3, 'ko', ...
     'MarkerFaceColor', 'k')
grid on
xlabel('Frequency (Hz)')
ylabel('Peak Amplitude (mm)')
subplot(2,1,2)
for i=1:length(ns)
    plot(t, uouts{i}(:, mi), 'LineWidth', 2); hold on
end
plot(xlim, gap*[1 1], 'k--'); hold on
plot(xlim, -gap*[1 1], 'k--'); hold on
grid on
xlabel('Scaled Time')
ylabel('Response')
drawnow

%% Animate Deflection
sc = 1;
U = UCs{fi}(1:end-1, mi);
% U(MDL.Ndofs+(1:4*MDL.Ndofs)) = 0;  % Removing the first harmonic
% sc = 4;
U = Lb(1:2:end,:)*AFT(reshape(U, MDL.Ndofs, Nhc)', h,Nt, 'f2t')';
if animfig
    figure(101);
    clf()
    ncyc = 1;
    for ti=repmat(1:16:Nt, 1,ncyc)
        clf()
        subplot(2,1,1)
        plot(Xn1(:,1), Xn1(:,2), 'k--'); hold on
        plot(Xn2(:,1), Xn2(:,2), 'k--');
        plot(Xn3(:,1), Xn3(:,2), 'k--');
        
        plot(Xn1(:,1), Xn1(:,2)+sc*U(1:Nn1,ti), 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'w');
        plot(Xn2(:,1), Xn2(:,2)+sc*U(Nn1+(1:Nn2),ti), 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'w');
        plot(Xn3(:,1), Xn3(:,2)+sc*U(Nn1+Nn2+(1:Nn2),ti), 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'w');

        ps = [Xn2(1,2)+sc*U(Nn1+1,ti) Xn1(Nn1,2)+sc*U(Nn1,ti)];
        if diff(ps)>0
            plot(Xn1(Nn1, 1)*[1 1], ps, 'ro-', 'LineWidth', 2, ...
                 'MarkerFaceColor', 'r');
        end
        ps = [Xn1(Nn1,2)+sc*U(Nn1,ti) Xn3(1,2)+sc*U(Nn1+Nn2+1,ti)];
        if diff(ps)>0
            plot(Xn1(Nn1, 1)*[1 1], ps, 'ro-', 'LineWidth', 2, ...
                 'MarkerFaceColor', 'r');
        end    
        axis tight
        ylim(max(abs(uouts{1}(:, mi))+gap)*[-1 1])
        grid on

        nsp = [Nn1 Nn1+1 Nn1+Nn2+1];
        gps = [0 gap -gap];
        subplot(2,1,2)
        for i=1:length(ns)
            plot(t, gps(i)+sc*U(nsp(i),:), 'LineWidth', 2); hold on
            plot(t(ti), gps(i)+sc*U(nsp(i),ti), 'ko', 'MarkerFaceColor', 'k');
        end
        plot(xlim, gap*[1 1], 'k--'); hold on
        plot(xlim, -gap*[1 1], 'k--'); hold on
        grid on
        xlabel('Scaled Time')
        ylabel('Response')
        axis tight
        drawnow
    end
end

