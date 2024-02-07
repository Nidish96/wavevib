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
bcs = [struct('i', 1, 'cofs', [1 1 1 1; 1 -1 1j -1j], ...
              'cofs0', [eye(2) zeros(2)]);  % fixed end
       struct('i', 5, 'cofs', [1 1 1 1; 1 -1 1j -1j], ...
              'cofs0', [eye(2) zeros(2)]);     % fixed end
       struct('i', 7, 'cofs', [1 1 1 1; 1 -1 1j -1j], ...
              'cofs0', [eye(2) zeros(2)]);     % fixed end
       struct('i', 3, 'cofs', [1 1 -1 -1], ...
              'cofs0', [0 0 1 0]);        % Moment-free
       struct('i', 4, 'cofs', [1 1 -1 -1], ...
              'cofs0', [0 0 1 0]);        % Moment-free
       struct('i', 6, 'cofs', [1 1 -1 -1], ...
              'cofs0', [0 0 1 0])];       % Moment-free


%% Setup Excitation
Mx = @(w,xi) inv([Ey*Iy*Klib.K(w,xi)^3*[-1 1 1j -1j];
                  Ey*Iy*Klib.K(w,xi)^2*[-1 -1 1 1];
                  [1 1 1 1];
                  Klib.K(w,xi)*[1 -1 1j -1j]]);
excs = struct('i', 2, 'nh', 1, ...
              'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0], ...
              'rcofs0', [0;0;0;1/(6*Ey*Iy)]);
%'nh' sets the harmonic at which to apply the excitation

%% Setup AFT parameters
h = (0:16)'; % (0:16)'
Nt = 2^10;

%% Setup Joints
cofs = [-[zeros(1,4) 1 -1 -1j 1j zeros(1,4)];
        -[zeros(1,8) 1 -1 -1j 1j];
        kron([1 -1 -1], [1 -1 -1j 1j])];
cofs0 = [-[zeros(1,4) 0 0 0 6 zeros(1,4)];
         -[zeros(1,8) 0 0 0 6];
         kron([1 -1 -1], [0 0 0 1])];
joint_nl = struct('type', 3, 'is', [3 4 6], ...
                  'cofs', @(w,xi) (Ey*Iy*Klib.K(w,xi)^3)*cofs(1:2,:), ...
                  'nl', @(Uw) HCONTACT(Uw, knl, gap, h, Nt), ...
                  'nldcofs', kron([1 -1 0;-1 0 1], [1 1 1 1]), ...
                  'nlfcofs', [1 0;0 -1], ...
                  ...
                  'cofs0', (Ey*Iy)*cofs0(1:2,:), ...
                  'nldcofs0', kron([1 -1 0;-1 0 1], [1 0 0 0]), ...
                  'nlfcofs0', [1 0;0 -1]);
joint_l = struct('type', 3, 'is', [3 4 6], 'cofs', cofs(3,:), 'cofs0', cofs0(3,:));

joints(2) = joint_l;
fnms = fieldnames(joint_nl);
for i=1:length(fnms)
    joints(1).(fnms{i}) = joint_nl.(fnms{i});
    % joints(1) = setfield(joints(1), fnms{i}, getfield(joint_nl, fnms{i}));
end

% joints(2) = joint_nl;
% fnms = fieldnames(joint_l);
% for i=1:length(fnms)
%     joints(1) = setfield(joints(1), fnms{i}, getfield(joint_l, fnms{i}));
% end


% cofs = @(w,xi) (Ey*Iy*Klib.K(w,xi)^3)*[-[zeros(1,4) 1 -1 -1j 1j zeros(1,4)];
%                                        -[zeros(1,8) 1 -1 -1j 1j];
%                                        kron([1 -1 -1], [1 -1 -1j 1j])];
% cofs0 = @(w,xi) (Ey*Iy)*[-[zeros(1,4) 0 0 0 6 zeros(1,4)];
%                          -[zeros(1,8) 0 0 0 6];
%                          kron([1 -1 -1], [0 0 0 1])];
% joints = struct('type', 3, 'is', [3 4 6], ...
%                 'cofs', cofs, 'nl', @(Uw) HCONTACT(Uw, knl, gap, h, Nt), ...
%                 'nldcofs', kron([1 -1 0;-1 0 1], [1 1 1 1]), ...
%                 'nlfcofs', [1 0;0 -1;0 0], ...
%                 ...
%                 'cofs0', cofs0, ...
%                 'nldcofs0', kron([1 -1 0;-1 0 1], [1 0 0 0]), ...
%                 'nlfcofs0', [1 0;0 -1;0 0]);

% cofs = [-[zeros(1,4) 1 -1 -1j 1j zeros(1,4)];
%         -[zeros(1,8) 1 -1 -1j 1j];
%         kron([1 -1 -1], [1 -1 -1j 1j])];
% cofs0 = [-[zeros(1,4) 0 0 0 6 zeros(1,4)];
%          -[zeros(1,8) 0 0 0 6];
%          kron([1 -1 -1], [0 0 0 1])];
% joints = struct('type', 3, 'is', [3 4 6], ...
%                 'cofs', cofs, 'nl', @(Uw) HCONTACT(Uw, knl, gap, h, Nt), ...
%                 'nldcofs', kron([1 -1 0;-1 0 1], [1 1 1 1]), ...
%                 'nlfcofs', @(w,xi) [1 0;0 -1;0 0]/(Ey*Iy*Klib.K(w,xi)^3), ...
%                 ...
%                 'cofs0', cofs0, ...
%                 'nldcofs0', kron([1 -1 0;-1 0 1], [1 0 0 0]), ...
%                 'nlfcofs0', [1 0;0 -1;0 0]/(Ey*Iy));

%% Pre-Processing
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);
% Klib.dKdw = @(w,xi) Klib.dKdw(w+eps,xi);  % Adjusting for w=0.

%% Conduct Nonlinear Forced Response Analysis using Continuation
Npts = pcs(end).irange(end);
Nwc = size(wcomps,1);

Nh = length(h);
Nhc = sum((h==0)+2*(h~=0));
[zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);

Wst = 2*pi*4;
Wen = 2*pi*12;

Famps = [0.4, 1.0 2.5, 5.0];

fi = 3;

[Amat, dAmatdw, ~, Fv, ~, ~, JEV] = WVAMAT([Wst;0], h, pcs, bcs, joints, Klib, 'r');
ari0_o = Amat\(Fv*Famps(fi));

opt = optimoptions('fsolve', 'specifyObjectiveGradient', true, 'Display', 'iter');
tic
aris_o = fsolve(@(ari) WVHBRESFUN([ari;Wst], Famps(fi), h, pcs, bcs, joints, Klib), ...
                ari0_o, opt);
toc

%% Test WVAMATr
[Amatr, dAmatdwr, dAmatdxir, Fvr, dFvdwr, dFvdxir, JEVr, RECOV] = WVAMATr([Wst;0], h, pcs, bcs, joints, Klib, 'r');

ari0 = Amatr\(Fvr*Famps(fi));
tic
aris = fsolve(@(ari) WVHBRESFUNr([ari;Wst], Famps(fi), h, pcs, bcs, joints, Klib), ...
              ari0, opt);
toc

norm(aris_o)
norm(aris_o-(RECOV.Qm*aris+RECOV.Qv*Famps(fi)))

nrep = 4;
tic
for i=1:nrep
    WVHBRESFUNr([aris;Wst], Famps(fi), h, pcs, bcs, joints, Klib);
end
toc/nrep

% %% Continuation: Full Representation
% Copt = struct('Nmax', 200, 'angopt', 5e-2, 'DynDscale', 1, 'solverchoice', 2);
% [Amat, ~, ~, Fv] = WVAMAT([Wst;0], h, pcs, bcs, joints, Klib, 'r');

% ari0 = Amat\(Fv*Famps(fi));
% Copt.Dscale = [1e-3*ones(size(ari0));Wst];

% dw = 0.75;
% tic
% ariwC = CONTINUE(@(ariw) WVHBRESFUN(ariw, Famps(fi), h, pcs, bcs, joints, Klib), ...
%                  ari0, Wst, Wen, dw, Copt);
% toc

% %% Continuation: Reduced Representation
% Copt = struct('Nmax', 20, 'angopt', 5e-2, 'DynDscale', 1, 'solverchoice', 2);
% [Amatr, ~, ~, Fvr] = WVAMATr([Wst;0], h, pcs, bcs, joints, Klib, 'r');

% ari0 = Amatr\(Fvr*Famps(fi));
% Copt.Dscale = [1e-3*ones(size(ari0));Wst];

% dw = 0.75;
% tic
% ariwC = CONTINUE(@(ariw) WVHBRESFUNr(ariw, Famps(fi), h, pcs, bcs, joints, Klib), ...
%                  ari0, Wst, Wen, dw, Copt);
% toc

%% Continuation: Reduced Representation
hsc = kron(h, [1;1]);
if h(1)==0
    hsc = hsc(2:end);
end

Copt = struct('Nmax', 150, 'angopt', 5e-2, 'DynDscale', 1, 'solverchoice', 2);
[Amatr, ~, ~, Fvr, ~, ~, JEVr] = WVAMATr([Wst;0], h, pcs, bcs, joints, Klib, 'r');

ari0 = Amatr\(Fvr*Famps(fi));
% Copt.Dscale = [1e0*kron(1./max(hsc,1), ones(2,1));Wst];
% Copt.Dscale = [1e-2./sum(abs(Amatr),2); Wst];
Copt.Dscale = [ones(size(ari0))*1e-2; Wst];

dw = 0.005;  % 0.01
tic
ariwC = CONTINUE(@(ariw) WVHBRESFUNr(ariw, Famps(fi), h, pcs, bcs, joints, Klib), ...
                 ari0, Wst, Wen, dw, Copt);
toc

% dw = 0.25;
% Sopt = struct('stepmax', 100, 'jac', 'x', 'dynamicDscale', 1);
% Sopt.Dscale = [ones(size(ari0))*1e-2; Wst];
% tic
% ariwC = solve_and_continue(ari0, @(ariw) WVHBRESFUNr(ariw, Famps(fi), h, pcs, bcs, ...
%                                                      joints, Klib), ...
%                            Wst, Wen, dw, Sopt);
% toc

% Stability
pds = zeros(size(ariwC,2), 1);

AriwC = zeros(Npts*Nwc*Nhc+1, size(ariwC,2));
for i=1:size(ariwC,2)
    W = ariwC(end,i);

    [pAmat, RECOV] = WBPERJACFUNr(W, ariwC(:, i), Famps(fi), h, pcs, bcs, joints, Klib);
    pds(i) = det(pAmat);
    AriwC(:, i) = [RECOV.Qm*ariwC(1:end-1,i)+RECOV.Qv*Famps(fi); ariwC(end,i)];
end

%%
stabi = sign(conv(sign(pds(:)')==sign(pds(1)), [1 1], 'same'));
opi = (13:16)';

Rout0 = ((1:Npts*Nwc)==opi(1));
Rout = 2*sum((1:Npts*Nwc)==opi);

Rv = blkdiag(Rout0, kron(speye(Nhc-1), Rout));
Aout = Rv*AriwC(1:end-1,:);

figure(1)
% clf()
plot(AriwC(end,:)./stabi/2/pi, abs([1 1j]*Aout(2:3,:)), '.-'); hold on
plot(AriwC(end,:)./(1-stabi)/2/pi, abs([1 1j]*Aout(2:3,:)), '.--')
xlim([Wst Wen]/2/pi)

xlabel('Frequency (Hz)')
ylabel('Response (m)')
