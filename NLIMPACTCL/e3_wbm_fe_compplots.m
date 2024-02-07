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

savfig = true;
animfig = false;

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
cofs0 = @(xi) [-[zeros(1,4) 0 0 0 6 zeros(1,4)];
               -[zeros(1,8) 0 0 0 6];
               kron([1 -1 -1], [0 0 0 1])];
joints = struct('type', 3, 'is', [3 4 6], ...
                'cofs', cofs, 'nl', [], ...  % Setup nl after settin 'h'
                'nldcofs', @(w,xi) kron([1 -1 0;-1 0 1], [1 1 1 1]), ...
                'nlfcofs', @(w,xi) [1 0;0 -1;0 0]/(Ey*Iy*Klib.K(w,xi)^3), ...
                ...
                'cofs0', cofs0, ...
                'nldcofs0', @(xi) kron([1 -1 0;-1 0 1], [1 0 0 0]), ...
                'nlfcofs0', @(xi) [1 0;0 -1;0 0]/(Ey*Iy));

%% Pre-Processing
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);

%% Conduct Nonlinear Forced Response Analysis
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

h = (0:16)';  % (0:16)'
Nt = 2^10;
Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);
[zinds1,hinds1,rinds01,rinds1,iinds1] = HINDS(1, h);

load(sprintf('./DATS/E_NLIMPDAT_2Knl%d_H%d.mat', knl,Nh), ...
     'acC', 'pds', 'h', 'Famps', 'Klib', 'wcomps', ...
     'pcs', 'bcs', 'joints', 'excs', ...
     'knl', 'gap', 'Wst', 'Wen')

Nf = length(Famps);
joints.nl = @(Uw) HCONTACT(Uw, knl, gap, h, Nt);

stabis = cellfun(@(pp) sign(conv(sign(pp(:)')==sign(pp(1)), [1 1], 'same')), ...
                 pds, 'Uniformoutput', false);

%% Load FE Results
Nn1 = 10;
hfe = h;

% Nn1 = 20;
% hfe = (0:3)';

fe = load(sprintf('./DATS/E1_FERES_2Knl%d_N%d_H%d.mat', 2*knl, Nn1, max(hfe)), ...
          'UCs', 'Fl', 'Famps', 'Rb', 'h', 'Nt', 'Wst', 'Wen', 'Nn1', 'Lb');
fe.Rrel = fe.Lb((fe.Nn1-1)*2+1,:)-fe.Lb((fe.Nn1+1-1)*2+1,:);
fe.Nh = length(fe.h);
fe.Nhc = sum((fe.h==0)+2*(fe.h~=0));
fe.Lhc = zeros(Nh, fe.Nhc);
fe.Lhc(1,1) = 1;
fe.Lhc(2:fe.Nh, [2:2:fe.Nhc 3:2:fe.Nhc]) = kron([1 1j], eye(fe.Nh-1));

%% Plot Forced Response
opi1 = (13:16)';  % x1
opi2 = (17:20)';  % x2
opi3 = (25:28)';  % x3
                  % hi = find(h==1);  % Choose which harmonic to plot
hsel = (h==1)';
% hsel = 2*[1 0.5*h(h~=0)'];

% fe.Rout = fe.Rb;
% wbm.Rout = ((1:Npts*Nwc)==opi1(1));
% wbm.Rout = 2*sum((1:Npts*Nwc)==opi1);  % Abs response at x1

fe.Rout = fe.Rrel;
wbm.Rout0 = ((1:Npts*Nwc)==opi1(1))-((1:Npts*Nwc)==opi2(1));
wbm.Rout = 2*(sum((1:Npts*Nwc)==opi1)-sum((1:Npts*Nwc)==opi2));  % relative coords

colos = DISTINGUISHABLE_COLORS(length(Famps))*0.75;
npfe = 20;

for hi=[1]   
    hsel = (h==hi)';
    figure(hi+1)
    set(gcf, 'Color', 'white')
    clf()
    aa = gobjects(size(fe.Famps));
    for fi=1:length(fe.Famps)
        if fi<=length(fe.Famps)
            Uh = kron(fe.Lhc, fe.Rout)*fe.UCs{fi}(1:end-1,:);
        end
        aout = blkdiag(wbm.Rout0, kron(eye(Nh-1), wbm.Rout))*acC{fi}(1:end-1,:);
        
        % subplot(2,1,1)
        aa(fi)=plot(acC{fi}(end,:)./stabis{fi}/2/pi, abs(hsel*aout)*1e3, ...
                    '-', 'LineWidth', 2, 'Color', colos(fi,:)); hold on
        legend(aa(fi), sprintf('F = %.1f N', Famps(fi)));
        plot(acC{fi}(end,:)./(1-stabis{fi})/2/pi, abs(hsel*aout)*1e3, ...
             '--', 'LineWidth', 2, 'Color', colos(fi,:))
        if fi<=length(fe.Famps)
            phs = angle(Uh(h==1,:));
            minds = fix(interp1(phs, 1:size(Uh,2), linspace(phs(2), phs(end-1), npfe)));
            
            bb = plot(fe.UCs{fi}(end,:)/2/pi, abs(hsel*Uh*1e3), ...
                      'o-', 'LineWidth', 0.25, 'Color', colos(fi,:), ...
                      'MarkerFaceColor', 'w', 'MarkerSize', 4, ...
                      'MarkerIndices', minds);
            legend(bb, 'FE-HB')
        end
        % subplot(2,1,2)
        % plot(acC{fi}(end,:)/2/pi, rad2deg(angle(aout)), ...
        %      '-', 'LineWidth', 1.5); hold on
        % if fi<=length(fe.Famps)
        %     plot(fe.UCs{fi}(end,:)/2/pi, rad2deg(angle([1 1j]*Uh(2:3,:))), ...
        %          '--', 'LineWidth', 2, 'Color', colos(fi,:))
        % end
        % grid on
        % ylabel('Phase (degs)')
    end
    grid on
    xlabel('Excitation Frequency (Hz)')
    ylabel(sprintf('H%d Response (mm)', hi))
    % subplot(2,1,1)
    if hsel(h==1)
        plot(xlim, gap*[1 1]*1e3, 'k-.')
        text(10, gap*1.1*1e3, 'Gap limit')
    end
    xlim(sort([Wst Wen]/2/pi))
    ll = legend([aa bb], 'Location', 'northwest');
    drawnow

    if hi==1
        % Customize 4th entry
        tm = ll.EntryContainer.NodeChildren(1);
        iconset = tm.Icon.Transform.Children.Children;
        iconset(1).VertexData(1) = 0.85;
        iconset(1).LineWidth = 1;
        iconset(2).ColorData = uint8([0 0 0 255]');

        ic1 = copy(iconset);
        ic1(1).Parent = iconset(1).Parent;
        ic1(1).EdgeColorData = uint8([colos(1,:) 1]'*255);
        ic1(1).LineWidth = 1;
        ic1(1).VertexData(1) = 0.15;

        ic2 = copy(iconset);
        ic2(1).Parent = iconset(1).Parent;
        ic2(1).EdgeColorData = uint8([colos(2,:) 1]'*255);
        ic2(1).VertexData(1) = 0.5;
        ic2(1).LineWidth = 1;

        ll = legend(aa, 'Location', 'northwest');
    else
        ll.Visible = 'off';
    end
    % drawnow
    if hi==0
        xlim([7.5 10.1])
    end
    
    if savfig
        savefig(sprintf('./FIGS/E3_FECOMP_H%d.fig', hi))
        export_fig(sprintf('./FIGS/E3_FECOMP_H%d.png', hi), '-dpng', '-r300')
        export_fig(sprintf('./FIGS/E3_FECOMP_H%d.eps', hi), '-depsc')       
    end
end

% subplot(2,1,2)
% xlim(sort([Wst Wen]/2/pi))
% set(gca, 'YTick', -180:45:180)
% ylim([0 1]*180)
% xlabel('Frequency (Hz)')

%% Plot RMS Response
hsel = [2 ones(1, Nh-1)];

figure(3)
set(gcf, 'Color', 'white')
clf()
for pli=1:2
    aa = gobjects(size(fe.Famps));
    for fi=1:length(fe.Famps)
        if fi<=length(fe.Famps)
            Uh = kron(fe.Lhc, fe.Rout)*fe.UCs{fi}(1:end-1,:);
        end
        aout = blkdiag(wbm.Rout0, kron(eye(Nh-1), wbm.Rout))*acC{fi}(1:end-1,:);
        
        % subplot(2,1,1)
        aa(fi)=plot(acC{fi}(end,:)./stabis{fi}/2/pi, abs(hsel*aout)*1e3, ...
                    '-', 'LineWidth', 2, 'Color', colos(fi,:)); hold on
        legend(aa(fi), sprintf('F = %.1f N', Famps(fi)));
        plot(acC{fi}(end,:)./(1-stabis{fi})/2/pi, abs(hsel*aout)*1e3, ...
             '--', 'LineWidth', 2, 'Color', colos(fi,:))
        if fi<=length(fe.Famps)
            phs = angle(Uh(h==1,:));
            minds = fix(interp1(phs, 1:size(Uh,2), linspace(phs(2), phs(end-1), npfe)));
            
            bb = plot(fe.UCs{fi}(end,:)/2/pi, abs(hsel*Uh*1e3), ...
                      'o-', 'LineWidth', 2, 'Color', colos(fi,:), ...
                      'MarkerFaceColor', 'w', 'MarkerSize', 4, ...
                      'MarkerIndices', minds, 'LineWidth', 0.25);
            legend(bb, 'FE-HB')
        end
    end
    grid on
    xlabel('Excitation Frequency (Hz)')
    ylabel('$\sqrt{2}\times$RMS Response (mm)')
    if hsel(h==1)
        plot(xlim, gap*[1 1]*1e3, 'k-.')
        tx=text(10, gap*1.1*1e3, 'Gap limit');
    end
    xlim(sort([Wst Wen]/2/pi))
    ll = legend([aa bb], 'Location', 'northwest');
    ll.Visible = 'off';

    xll = [9.1 9.42];
    yll = [3.4 3.9];
    if pli==1
        plot(xll([1 2 2 1 1]), yll([1 1 2 2 1]), 'k-')
        annotation('arrow', [.62 .45], [.57 .67]);
        
        axes('Position', [.2 .65 .25 .25]);
    else
        xlim(xll);
        ylim(yll)
        xlabel(''); ylabel('');
        tx.Visible='off';
    end
end

if savfig
    savefig('./FIGS/E3_FECOMP_RMSA.fig')
    export_fig('./FIGS/E3_FECOMP_RMSA.eps', '-depsc');
end

%% Plot of Higher harmonics

% his = [2 3 5 7 8 11];
his = [2 3 4 5 7 8];

figure(4)
poss = get(gcf, 'Position');
set(gcf, 'Color', 'white', 'Position', [poss(1:2) [720 450]*1.5]);
clf()
tl=tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% tl.TileSpacing = 'compact';
for hi=1:length(his)
    hsel = (h==his(hi))';
    % subplot(2, 3, hi)
    nexttile;
    aa = gobjects(size(fe.Famps));
    for fi=1:length(fe.Famps)
        if fi<=length(fe.Famps)
            Uh = kron(fe.Lhc, fe.Rout)*fe.UCs{fi}(1:end-1,:);
        end
        aout = blkdiag(wbm.Rout0, kron(eye(Nh-1), wbm.Rout))*acC{fi}(1:end-1,:);
        
        % subplot(2,1,1)
        aa(fi)=plot(acC{fi}(end,:)./stabis{fi}/2/pi, abs(hsel*aout)*1e3, ...
                    '-', 'LineWidth', 2, 'Color', colos(fi,:)); hold on
        legend(aa(fi), sprintf('F = %.1f N', Famps(fi)));
        plot(acC{fi}(end,:)./(1-stabis{fi})/2/pi, abs(hsel*aout)*1e3, ...
             '--', 'LineWidth', 2, 'Color', colos(fi,:))
        if fi<=length(fe.Famps)
            phs = angle(Uh(h==1,:));
            minds = fix(interp1(phs, 1:size(Uh,2), linspace(phs(2), phs(end-1), npfe)));
            
            bb = plot(fe.UCs{fi}(end,:)/2/pi, abs(hsel*Uh*1e3), ...
                      'o-', 'LineWidth', 2, 'Color', colos(fi,:), ...
                      'MarkerFaceColor', 'w', 'MarkerSize', 4, ...
                      'MarkerIndices', minds, 'LineWidth', 0.25);
            legend(bb, 'FE-HB')
        end
    end
    grid on
    if hi>3
        xlabel('Frequency (Hz)')
    end
    ttl=title(sprintf('\\textbf{Harmonic %d}', his(hi)));
    % ttl.Position(2) = ttl.Position(2)*0.9;
    if hi==1 || hi==4
        ylabel('Response (mm)')
    end
    % subplot(2,1,1)
    xlim([7.5 10.1])
    ll = legend([aa bb], 'Location', 'northwest');
    ll.Visible = 'off';
end

if savfig
    savefig('./FIGS/E3_FECOMP_HARMS.fig');
    export_fig('./FIGS/E3_FECOMP_HARMS.eps', '-depsc');
end


% %% Plot the response at one point (From WBM Solution)
% t = linspace(0, 2*pi, Nt+1)'; t=t(1:end-1);
% fi = 2;
% opis = {(13:16), (17:20), (25:28)};  % Points x1, x2, x3 in paper

% Uouts = cell(size(opis));
% uouts = cell(size(opis));
% for i=1:length(opis)
%     Rout = speye(Npts*Nwc);
%     Rout = kron(eye(Nh), 2*sum(Rout(opis{i},:)));
%     if h(1)==0
%         Rout(1, opis{i}) = [1 0 0 0];
%     end
    
%     Uouts{i} = Rout*acC{fi}(1:end-1,:);
%     Uouth = zeros(Nhc, size(Uouts{i},2));
%     Uouth([rinds01 rinds1 iinds1], :) = ...
%         [real(Uouts{i}(zinds1,:));
%          real(Uouts{i}(hinds1,:)); imag(Uouts{i}(hinds1,:))];
%     uouts{i} = AFT(Uouth, h,Nt, 'f2t');
%     if i>1
%         uouts{i} = uouts{i}+(-1)^i*gap;
%     end
% end
% [~, mi] = max(max(abs(uouts{1})));
% % [~, mi] = min(abs(acC{fi}(end,:)/2/pi-9));

% colos = DISTINGUISHABLE_COLORS(length(opis))*0.75;
% aa = gobjects(size(opis));
% figure(2)
% set(gcf, 'Color', 'white')
% clf()
% subplot(2,1,1)
% for i=1:length(opis)
%     aa(i) = plot(acC{fi}(end,:)/2/pi, max(abs(uouts{i}))*1e3, '-', ...
%                  'LineWidth', 1.5, 'Color', colos(i,:)); hold on
%     plot(acC{fi}(end,mi)/2/pi, max(abs(uouts{i}(:,mi)))*1e3, 'ko', ...
%          'MarkerFaceColor', 'k');

%     legend(aa(i), ['Point x' num2str(i)])
% end
% legend(aa, 'Location', 'northwest')
% axis tight
% grid on
% xlabel('Frequency (Hz)')
% ylabel('Peak Amplitude (mm)')
% subplot(2,1,2)
% for i=1:length(opis)
%     plot(t, uouts{i}(:, mi)*1e3, 'LineWidth', 2, ...
%          'Color', colos(i,:)); hold on
% end
% set(gca, 'XTick', 0:pi/2:2*pi, 'XTicklabel', ...
%          {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'})
% axis tight
% xlim([0 2*pi])
% plot(xlim, gap*[1 1]*1e3, 'k--'); hold on
% plot(xlim, -gap*[1 1]*1e3, 'k--'); hold on
% grid on
% xlabel('Scaled Time')
% ylabel('Response')
