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

Nn1 = 10;
Nn2 = 8;
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
h = (0:5)';
Nt = 128;
Nhc = sum((h==0)+2*(h~=0));
[~,~,zinds,rinds,iinds] = HINDS(1,h);

Fl = zeros(Nhc,1);
Fl(rinds(1)) = 1;
Fl = kron(Fl, Fv);

Famps = [0.4 1.0 2.5];  % 0.336, 0.336

Wst = 2*pi*4;
Wen = 2*pi*12;
dw  = 0.75;

Copt = struct('Nmax', 600, 'angopt', 2.5e-2, 'DynDscale', 1, ...
    'Dscale', [1e-6*ones(MDL.Ndofs*Nhc,1);Wst]);

UCs = cell(size(Famps));
for fi=1:length(Famps)
    U0 = HARMONICSTIFFNESS(MDL.M,MDL.C,MDL.K,Wst,h)\(Fl*Famps(fi));

    UCs{fi} = CONTINUE(@(Uw) MDL.HBRESFUN(Uw, Fl*Famps(fi), h, Nt), U0, Wst, Wen, dw, Copt);
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