clc
clear all
addpath('../ROUTINES/SOLVERS/')

%%
h = 1;
Nt = 128;

Ey = 190e9;
G = 77.5e9;
nu = 0.29;
rho = 7680;
wid = 0.1;
brd = 0.4;
Ar = wid*brd;
Iy = wid^3*brd/12;
L0 = 2.0;
kap = 10*(1+nu)/(12+11*nu);

% Dispersion Relationship
a = Ey*kap*G*Iy;
b = @(w, xi) -(w.^2*rho*Ar+rho*Ey*Iy);
c = @(w, xi) w.^4*rho^2*Ar-w.^2*rho*kap*G*Ar;
Klib = [struct('K', @(w,xi) sqrt((-b(w,xi)+sqrt(b(w,xi).^2-4*a*c(w,xi)))/(2*a))); ...
    struct('K', @(w,xi) sqrt((-b(w,xi)-sqrt(b(w,xi).^2-4*a*c(w,xi)))/(2*a)))];
wcomps = [1j 1;  % +jK1
         -1j 1;  % -jK1
          1j 2;  % +jK2
         -1j 2]; % -jK2

% The following represents the ratio of "theta"'s wave coefficient wrt
% "v"'s, i.e., the ratio of rotation wave coefficient w.r.t. translational
% wave coefficient
rel = struct('P', @(w,xi,cc,K) ((cc*K.K(w,xi)).^2*kap*G*Ar-w.^2*rho*Ar)./(-1j*cc*K.K(w,xi)*kap*G*Ar));
Ps = arrayfun(@(i) @(w,xi) rel.P(w,xi,wcomps(i,1),Klib(wcomps(i,2))), 1:4, 'UniformOutput', false);
Ks = arrayfun(@(i) Klib(wcomps(i,2)).K, 1:4, 'UniformOutput', false);

%%
% Pieces
pcs = [struct('coords', [0; L0/3; L0], 'wcomps', wcomps);
    struct('coords', [L0; 2*L0], 'wcomps', wcomps)];
% fix-fix
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1]);
    struct('i', 1, 'cofs', @(w,xi) [Ps{1}(w,xi) Ps{2}(w,xi) Ps{3}(w,xi) Ps{4}(w,xi)]);
    struct('i', 5, 'cofs', @(w,xi) [1 1 1 1]);
    struct('i', 5, 'cofs', @(w,xi) [Ps{1}(w,xi) Ps{2}(w,xi) Ps{3}(w,xi) Ps{4}(w,xi)])];
% % fix-free
% bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1]);
%     struct('i', 1, 'cofs', @(w,xi) [Ps{1}(w,xi) Ps{2}(w,xi) Ps{3}(w,xi) Ps{4}(w,xi)]);
%     struct('i', 5, 'cofs', @(w,xi) [1 1 1 1]);
%     struct('i', 5, 'cofs', @(w,xi) [Ps{1}(w,xi) Ps{2}(w,xi) Ps{3}(w,xi) Ps{4}(w,xi)])];


%% Excitation
Mx = @(w,xi) inv([1j*Ks{1}(w,xi)-Ps{1}(w,xi) 1j*Ks{2}(w,xi)-Ps{2}(w,xi) 1j*Ks{3}(w,xi)-Ps{3}(w,xi) 1j*Ks{4}(w,xi)-Ps{4}(w,xi);
    1j*Ks{1}(w,xi)*Ps{1}(w,xi) 1j*Ks{2}(w,xi)*Ps{2}(w,xi) 1j*Ks{3}(w,xi)*Ps{3}(w,xi) 1j*Ks{4}(w,xi)*Ps{4}(w,xi);
    1 1 1 1;
    Ps{1}(w,xi) Ps{2}(w,xi) Ps{3}(w,xi) Ps{4}(w,xi)]);

excs = struct('i', 2, 'nh', 1, ...
    'rcofs', @(w,xi) Mx(w,xi)*[1/(2*kap*G*Ar);0;0;0]); 

%% Joints
Vrow = @(w,xi) kap*G*Ar*[1j*Ks{1}(w,xi)-Ps{1}(w,xi) 1j*Ks{2}(w,xi)-Ps{2}(w,xi) 1j*Ks{3}(w,xi)-Ps{3}(w,xi) 1j*Ks{4}(w,xi)-Ps{4}(w,xi)];
Mrow = @(w,xi) - Ey*Iy*[1j*Ks{1}(w,xi)*Ps{1}(w,xi) 1j*Ks{2}(w,xi)*Ps{2}(w,xi) 1j*Ks{3}(w,xi)*Ps{3}(w,xi) 1j*Ks{4}(w,xi)*Ps{4}(w,xi)];
cofs = @(w,xi) [Vrow(w,xi) 0 0 0 0; Mrow(w,xi) 0 0 0 0;Vrow(w,xi) -Vrow(w,xi);Mrow(w,xi) -Mrow(w,xi)];
nldcofs = @(w,xi) [1 1 1 1 -1 -1 -1 -1; 
    kron([1 -1], [Ps{1}(w,xi) Ps{2}(w,xi) Ps{3}(w,xi) Ps{4}(w,xi)])];

kJs = diag([1e9 1e9]);
cJs = diag([320 420]);
gJs = diag([1e8 0]);
nljoints = struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs, ...
    'nl', @(Uw) HDUFF(Uw, kJs, cJs, gJs, h, Nt), ...
    'nlfcofs', @(w,xi) [eye(2); zeros(2)], ...
    'nldcofs', nldcofs);

% Fused Version
fusjoints = struct('type', 2, 'i', 3, 'j', 4);

% Linear Version
lincofs = @(w,xi) [(kJs(1,1)-cJs(1,1)*1j*w)*[1 1 1 1 -1 -1 -1 -1] + (kJs(1,2)-cJs(1,2)*1j*w)*kron([1 -1],[Ps{1}(w,xi) Ps{2}(w,xi) Ps{3}(w,xi) Ps{4}(w,xi)]);
    (kJs(2,1)-cJs(2,1)*1j*w)*[1 1 1 1 -1 -1 -1 -1] + (kJs(2,2)-cJs(2,2)*1j*w)*kron([1 -1],[Ps{1}(w,xi) Ps{2}(w,xi) Ps{3}(w,xi) Ps{4}(w,xi)])];
linjoints = struct( 'type', 2, 'i', 3, 'j', 4, ...
    'cofs', @(w,xi) cofs(w,xi)+[lincofs(w,xi);zeros(2,8)]);
%% Preprocessing
tic
[~, ~, fusjoints, ~, ~] = WBPREPROC(pcs, bcs, fusjoints, excs, Klib);
[~, ~, linjoints, ~, ~] = WBPREPROC(pcs, bcs, linjoints, excs, Klib);
[pcs, bcs, nljoints, excs, Klib] = WBPREPROC(pcs, bcs, nljoints, excs, Klib);
toc

%%
Nw = 1000;
Ws = linspace(eps,1e4,Nw);
Us = zeros(pcs(end).irange(end)*4, Nw);
Dv = zeros(size(Ws));
for iw=1:Nw
    [Amat, ~, ~, Fv] = WVAMAT([Ws(iw);0], h, pcs, bcs, linjoints, Klib);
    Us(:,iw) = Amat\Fv;
    Dv(iw) = det(Amat);
    fprintf('Done %d/%d\n', iw,Nw);
end

% [~, mi] = findpeaks(-abs(Dv));
% acms = zeros(Npts*Nwc*Nh,length(mi));
% wms = Ws(mi);
% opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', false, 'Display', 'iter');
% for i=1:length(mi)
%     Amat = WVAMAT([Ws(mi(i));0], h, pcs, bcs, linjoints, Klib);
% 
%     [U, S, V] = svd(Amat);
%     acms(:, i) = V(:,end);
% end
% 
% %% Plot mode shapes
% dL = [1 1 1 1];
% Nxp = 1000;
% Xs = linspace(0, 2*L0, Nxp);
% Vs = zeros(1,Nxp);
% k = 0;
% mi = 4;
% for n=1:length(pcs)
%     for i=1:size(pcs(n).coords,1)-1
%         k = k+1;
% 
%         inds = find(pcs(n).coords(i,1)<=Xs & pcs(n).coords(i+1,1)>=Xs);
%         
%         dXs = (Xs(inds)-pcs(n).coords(i,1));
% 
%         Vs(inds) = dL*(acms((k-1)*4+(1:4),mi).*exp([1j*Ks{1}(wms(mi),0)*dXs;-1j*Ks{1}(wms(mi),0)*dXs;1j*Ks{2}(wms(mi),0)*dXs;-1j*Ks{2}(wms(mi),0)*dXs]));
%     end
% end
% plot(Xs,Vs)

%%
Famp = 1e8*100;

Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));

Wst = 1150;
Wen = 1250;
dw = 0.1;

Copt = struct('Nmax', 300, 'angopt', 1e-1, 'DynDscale', 1);
[Amat, ~, ~, Fv] = WVAMAT([Wst;0], h, pcs, bcs, linjoints, Klib, 'r');
ari0 = Amat\Fv*Famp; 

ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famp, h, pcs, bcs, nljoints, Klib), ...
    ari0, Wst, Wen, dw, Copt);
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);
acC = zeros(Npts*Nwc*Nh+1, size(ariwC,2));
acC([zinds hinds end], :) = [ariwC(rinds0,:); ariwC(rinds,:)+1j*ariwC(iinds,:);ariwC(end,:)];
%% 
opi = 5:8;
figure(1)
% clf()
semilogy(Ws, abs(sum(2*Us(opi,:))), '.-'); hold on
aa=semilogy(acC(end,:), abs(sum(2*acC(opi,:)))/Famp, '.-'); hold on
% yyaxis right
% semilogy(Ws, abs(Dv))
% xlim([Wst Wen])
% xlim([200 203])
set(gca, 'YScale', 'linear')