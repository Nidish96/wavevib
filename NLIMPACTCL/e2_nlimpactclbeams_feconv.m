clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')
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
knl = 1210/2;    % 220,      880,    (242, 484, 880, 1210)
gap = 2.5e-3;   % 2.5e-3,   0.35,   0.35

Famps = [0.4 1.0 2.5];  % 0.336, 0.336
fi = 1;  % Chosen amplitude level (from above)

h = (0:16)';
Nt = 1024;
[rinds0,hinds,zinds,rinds,iinds] = HINDS(1,h);
Nhc = sum((h==0)+2*(h~=0));

%% Load Simulation with 10 elements (arbitrary number)
Nn1_ref = 10;
Nn2_ref = 8;

Xn1_ref = [linspace(0, L1, Nn1_ref)' zeros(Nn1_ref,1)];
Xn2_ref = [linspace(L1, L1+L2, Nn2_ref)' gap*ones(Nn2_ref,1)];
Xn3_ref = [linspace(L1, L1+L2, Nn2_ref)' -gap*ones(Nn2_ref,1)];

[~,mi]=min(abs(Xn1_ref(:,1)-xF1));
Xn1_ref(mi,1) = xF1;
Xn1_ref(:,1) = sort(Xn1_ref(:,1));

fe = load(sprintf('./DATS/E1_FERES_2Knl%d_N%d_H%d.mat', 2*knl, Nn1_ref, max(h)), ...
          'UCs', 'Fl', 'Famps', 'Rb', 'h', 'Nt', 'Wst', 'Wen', 'Nn1', 'Lb');

wbm = load('./DATS/E4_HCONV.mat', 'Hmaxs', 'Famps', 'aout', 'Wexcs');
Wexcs = wbm.Wexcs(1:3);

%% Compute solution at same frequency for a different mesh
Nn1s = [5 10 20 40 80 160];
% Nn1s = [5 10 20 40 80 120];
ttks = zeros(length(Nn1s), length(Famps));
Routs = zeros(Nhc, length(Nn1s), length(Famps));
Rrouts = zeros(Nhc, length(Nn1s), length(Famps));
opt = optimoptions('fsolve', 'specifyObjectiveGradient', true, 'Display', 'off');

if analyze
    for ni=1:length(Nn1s)
        Nn1 = Nn1s(ni);
        Nn2 = fix(.8*Nn1);
        Xn1 = [linspace(0, L1, Nn1)' zeros(Nn1,1)];
        Xn2 = [linspace(L1, L1+L2, Nn2)' gap*ones(Nn2,1)];
        Xn3 = [linspace(L1, L1+L2, Nn2)' -gap*ones(Nn2,1)];

        [~,mi]=min(abs(Xn1(:,1)-xF1));
        Xn1(mi,1) = xF1;
        Xn1(:,1) = sort(Xn1(:,1));
        iF1 = find(Xn1(:,1)==xF1);

        Xns = {Xn1; Xn2; Xn3};

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
        Nd = size(Lb,2);

        Mb = Lb'*Mf*Lb;
        Kb = Lb'*Kf*Lb;
        Cb = al*Mb + bt*Kb;  % Check this!
        Fv = Lb((iF1-1)*2+1,:)';
        Rb = Lb((Nn1-1)*2+1,:);  % Response pt -> last node of primary beam
        Rrel = Lb((Nn1-1)*2+1,:) - Lb((Nn1+1-1)*2+1,:); % Relative DOF

        %% Setup Nonlinearity
        nl1i = ([Nn1+1 Nn1]-1)*2+1;
        nl2i = ([Nn1 Nn1+Nn2+1]-1)*2+1;

        Lnls = [diff(Lb(nl1i,:));
                diff(Lb(nl2i,:))];

        nlfunc = @(t, u, ud) deal(max(knl*(u-gap),0), knl*(u>gap), zeros(size(u)));

        MDL = MDOFGEN(Mb,Kb,Cb,Lb);
        MDL = MDL.SETNLFUN(1+3, Lnls, nlfunc);
        clear Mb Kb Cb
        %% Setup HB
        Fl = zeros(Nhc,1);
        Fl(rinds(1)) = 1;
        Fl = kron(Fl, Fv);

        for fi=1:length(Famps)
            %% Interpolate initial condition
            Uh = kron(eye(Nhc), fe.Rb)*fe.UCs{fi}(1:end-1,:);
            Uh1 = abs([1 1j]*Uh(2:3,:));  % Peak first harmonic response
            [~, ipk] = max(Uh1);

            Uwpk_ref = fe.UCs{fi}(:, ipk);  % Solution at peak point
            Lb_ref = fe.Lb;
            Nd_ref = size(Lb_ref,2);

            UBms_ref = {Lb_ref(1:2*Nn1_ref,:)*reshape(Uwpk_ref(1:end-1), Nd_ref,Nhc);
                        Lb_ref(2*Nn1_ref+(1:2*Nn2_ref),:)*reshape(Uwpk_ref(1:end-1), Nd_ref,Nhc);
                        Lb_ref(2*Nn1_ref+2*Nn2_ref+(1:2*Nn2_ref),:)*reshape(Uwpk_ref(1:end-1), Nd_ref,Nhc)}; % {Bm1, Bm2, Bm3}
            Xns_ref = {Xn1_ref; Xn2_ref; Xn3_ref};
            W_ref = Uwpk_ref(end);

            UBms = {zeros(Nn1*2,Nhc); zeros(Nn2*2,Nhc); zeros(Nn2*2,Nhc)};
            for bi=1:3
                UBms{bi}(1:2:end,:) = interp1(Xns_ref{bi}(:,1), UBms_ref{bi}(1:2:end,:), ...
                                              Xns{bi}(:,1));
                UBms{bi}(2:2:end,:) = interp1(Xns_ref{bi}(:,1), UBms_ref{bi}(2:2:end,:), ...
                                              Xns{bi}(:,1));
            end

            U0 = reshape(Lb'*cell2mat(UBms), [], 1);  % Initial Guess for HB
            clear UBms

            % W_ref = Wexcs(fi);
	    Wexcs(fi) = W_ref;
            Usol = fsolve(@(U) MDL.HBRESFUN([U; W_ref], Fl*Famps(fi), h, Nt), U0, opt);

            nrep = 4;
            tic
            for nid=1:nrep
                MDL.HBRESFUN([Usol; W_ref], Fl*Famps(fi), h, Nt);
            end
            ttks(ni, fi) = toc/nrep;
            %%
            Routs(:, ni, fi) = (Rb*reshape(Usol,Nd,Nhc))';
            Rrouts(:, ni, fi) = (Rrel*reshape(Usol,Nd,Nhc))';
        end
        fprintf('Done %d/%d: %f\n', ni, length(Nn1s), ttks(ni, fi));
    end

    if savdat
        save('./DATS/E2_DATA.mat', 'Routs', 'Rrouts', 'Nn1s', 'Wexcs', 'ttks');
    end
else
    load('./DATS/E2_DATA.mat', 'Routs', 'Rrouts', 'Nn1s', 'Wexcs', 'ttks');
end

%% Load WBM Solutions
hci = find(wbm.Hmaxs==max(h));
wbmsol = squeeze(wbm.aout(1:(hci+1), hci, 1:3));
WBMSOL = zeros(Nhc, 3);
WBMSOL([zinds rinds iinds], :) = ...
    [wbmsol(rinds0, :); real(wbmsol(hinds, :)); imag(wbmsol(hinds, :))];

%% Plot Convergence
ndstot = Nn1s+2*fix(.8*Nn1s);

colos = DISTINGUISHABLE_COLORS(length(Famps))*0.75;

figure(1)
set(gcf, 'Color', 'white')
clf()
% errs = squeeze(rms(Rrouts-Rrouts(:,end,:))./rms(Rrouts(:,end,:)));
errs = squeeze(rms((Rrouts-Rrouts(:,end,:))./(Rrouts(:,end,:))));
% errs = [squeeze(rms(diff(Rrouts,[],2))./rms(Rrouts(:,end,:))); zeros(1,3)]; % Relative Cauchy Errors
% errs = squeeze(rms(Rrouts-permute(WBMSOL, [1 3 2])*2)./rms(2*permute(WBMSOL, [1 3 2])));
for fi=1:length(Famps)
    loglog(ndstot, errs(:, fi), 'o-', 'LineWidth', 2, ...
           'MarkerFaceColor', 'w', 'Color', colos(fi,:)); hold on
end
xl = xlim;
yl = ylim;
ndp = [1 1e3];
plot(ndp, 1e3./ndp.^4, 'k-.', 'LineWidth', 1);
xlim(xl)
ylim(yl)
grid on
grid minor
grid minor
legend([arrayfun(@(f) sprintf('F=%.2f N', f), Famps, 'UniformOutput', false), ...
        '$\mathcal{O}(N^{-4})$'], 'Location', 'northeast')
xlabel('Total No. of Nodes')
ylabel('Relative RMSE (all harmonics)')
if savfig
    savefig('./FIGS/E2_FECONV.fig')
    export_fig('./FIGS/E2_FECONV.eps', '-depsc');
end

figure(2)
set(gcf, 'Color', 'white')
clf()
for fi=1:length(Famps)
    loglog(ndstot, ttks(:, fi), 'o-', 'LineWidth', 2, ...
           'MarkerFaceColor', 'w', 'Color', colos(fi,:)); hold on
end
xl = xlim;
yl = ylim;
ndp = [1 1e3];
plot(ndp, ndp.^2/2e4, 'k-.', 'LineWidth', 1);
xlim(xl)
ylim(yl)
grid on
grid minor
grid minor
legend([arrayfun(@(f) sprintf('F=%.2f N', f), Famps, 'UniformOutput', false), ...
        '$\mathcal{O}(N^2)$'], 'Location', 'southeast')
xlabel('Total No. of Nodes')
ylabel('Avg. CPU Time (s)')

if savfig
    savefig('./FIGS/E2_FETTK.fig')
    export_fig('./FIGS/E2_FETTK.eps', '-depsc');
end
