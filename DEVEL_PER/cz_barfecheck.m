% clc
clear all
addpath('../ROUTINES/HARMONIC')
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/')

%%
Ey = 2.62e11;
rho = 1280;
Ar = 31.75e-3*54e-3;

xF = 0.28;
xJ = 1/3;
xL = 1.0;

kJ = 1e9;
cJ = 320;
gJ = 1e8;
fnl = @(t,u,ud) deal(gJ*u.^3,3*gJ*u.^2,zeros(size(u)));

Me = @(Le) rho*Ar*Le/6*[2 1;1 2];
Ke = @(Le) Ar*Ey/Le*[1 -1;-1 1];

Ne1 = 240;
Xn1 = linspace(0, xJ, Ne1+1);
[~,iF] = min(abs(Xn1-xF));
Xn1(iF) = xF; 

Ne2 = 2*Ne1;
Xn2 = linspace(xJ, xL, Ne2+1);

M1 = zeros(Ne1+1);
K1 = zeros(Ne1+1);
for ei=1:Ne1
    M1(ei:ei+1,ei:ei+1) = M1(ei:ei+1,ei:ei+1) + Me(diff(Xn1(ei:ei+1)));
    K1(ei:ei+1,ei:ei+1) = K1(ei:ei+1,ei:ei+1) + Ke(diff(Xn1(ei:ei+1)));
end
M2 = zeros(Ne2+1);
K2 = zeros(Ne2+1);
for ei=1:Ne2
    M2(ei:ei+1,ei:ei+1) = M2(ei:ei+1,ei:ei+1) + Me(diff(Xn2(ei:ei+1)));
    K2(ei:ei+1,ei:ei+1) = K2(ei:ei+1,ei:ei+1) + Ke(diff(Xn2(ei:ei+1)));
end

M = blkdiag(M1,M2);
K = blkdiag(K1,K2);
Lb = eye(Ne1+Ne2+2);
Lb = Lb(:,2:end-1);

Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;
Rb = Lb(Ne1+1+2,:);
LJ = Lb(Ne1+2,:)-Lb(Ne1+1,:);  % Relative Displacement at Joint
Fv = Lb(iF,:)';   % Forcing Vector

%% Setup MDOFGEN
GM = MDOFGEN(Mb,Kb+LJ'*kJ*LJ,LJ'*cJ*LJ,eye(size(Mb)));
GM = GM.SETNLFUN(1+3,LJ,fnl);

%% HB
Famp = 150e5;

Wst = 78e3;
Wen = 84e3;
dw = 0.1;

h = 1;
Nt = 128;

Nhc = sum((h==0)+2*(h~=0));
Fl = kron([zeros(h(1)==0,1); 1; 0;zeros(Nhc-2-(h(1)==0),1)], Fv);

Copt = struct('Nmax', 100, 'angopt', 2e-1, 'DynDscale', 1);
U0 = HARMONICSTIFFNESS(GM.M,GM.C,GM.K,Wst,h)\Fl*Famp;
UC = CONTINUE(@(Uw) GM.HBRESFUN(Uw,Fl*Famp,h,Nt), U0, Wst, Wen, dw, Copt);

%%
[zinds,hinds,rinds0,rinds,iinds] = HINDS(1, h);
Uoh = kron(eye(Nhc), Rb)*UC(1:end-1,:);
Uoc = zeros(length(h),size(Uoh,2));
Uoc([zinds hinds],:) = [Uoh(rinds0,:);Uoh(rinds,:)+1j*Uoh(iinds,:)];

figure(1)
% clf()
plot(UC(end,:)/1e3, abs(Uoc(h==1,:)), '.-'); hold on
xlim([80 81]);
ylim([0.05 0.3])