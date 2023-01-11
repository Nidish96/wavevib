clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',13)

%%
Ey = 190e9;
G = 77.5e9;
nu = 0.29;
rho = 7680;
wid = 0.2;
brd = 0.4;
Ar = wid*brd;
Iy = wid^3*brd/12;
L0 = 2.0;
kap = 10*(1+nu)/(12+11*nu);

Me = @(Le) rho*Ar*[(13*Le)/35, (11*Le^2)/210, (9*Le)/70, -(13*Le^2)/420;
    (11*Le^2)/210, Le^3/105, (13*Le^2)/420, -Le^3/140;
    (9*Le)/70, (13*Le^2)/420, (13*Le)/35, -(11*Le^2)/210;
    -(13*Le^2)/420, -Le^3/140, -(11*Le^2)/210, Le^3/105]; 
Ke = @(Le) Ey*Iy*[12/Le^3, 6/Le^2, -12/Le^3, 6/Le^2;
    6/Le^2, 4/Le, -6/Le^2, 2/Le;
    -12/Le^3, -6/Le^2, 12/Le^3, -6/Le^2;
    6/Le^2, 2/Le, -6/Le^2, 4/Le];

Ne1 = 20;
iJ = Ne1+1;
Xn1 = linspace(0, L0, Ne1+1);
[~, iF] = min(abs(Xn1-L0/3));
Xn1(iF) = L0/3;
M1 = zeros((Ne1+1)*2);
K1 = zeros((Ne1+1)*2);
for ei=1:Ne1
    is = (ei-1)*2+1;
    ie = (ei+1)*2;

    M1(is:ie, is:ie) = M1(is:ie, is:ie) + Me(diff(Xn1(ei:ei+1)));
    K1(is:ie, is:ie) = K1(is:ie, is:ie) + Ke(diff(Xn1(ei:ei+1)));
end

Ne2 = Ne1;
Xn2 = linspace(L0, 2*L0, Ne2+1);
M2 = zeros((Ne1+1)*2);
K2 = zeros((Ne1+1)*2);
for ei=1:Ne2
    is = (ei-1)*2+1;
    ie = (ei+1)*2;

    M2(is:ie, is:ie) = M2(is:ie, is:ie) + Me(diff(Xn2(ei:ei+1)));
    K2(is:ie, is:ie) = K2(is:ie, is:ie) + Ke(diff(Xn2(ei:ei+1)));
end

M = blkdiag(M1,M2);
K = blkdiag(K1,K2);
Lb = eye((Ne1+Ne2+2)*2);
Lb = Lb(:,3:end-2);
Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;
Fv = Lb((iF-1)*2+1,:)';
Rb = Fv';

%% Joint
kJs = diag([1e9 1e9]);
cJs = diag([320 320]);
gJ = 1e8;

LJ = Lb( (iJ-1)*2+(1:2), :)-Lb( (iJ)*2+(1:2), :);

KJ = LJ'*kJs*LJ;
CJ = LJ'*cJs*LJ;
fnl = @(t,u,ud) deal(gJ*u.^3, 3*gJ*u.^2, u*0);

%% 
GM = MDOFGEN(Mb, Kb+KJ, CJ, Lb);
GM = GM.SETNLFUN(1+3, LJ(1,:), fnl);

%%
[V,D] = eigs(GM.K, GM.M, 10, 'SM');
% [V,D] = eig(GM.K, GM.M);
Ws = sqrt(diag(D));
[~,iS] = sort(abs(Ws));
V = V(:,iS);
Ws = Ws(iS);

%%
Famp = 2e3*20; %[1, 10, 20]

h = 1;
Nt = 128;
Nhc = sum((h==0)+2*(h~=0));
Fl = kron([zeros(h(1)==0,1); 1; 0; zeros(Nhc-(h(1)==0)-2,1)], Fv);

Wst = 1055.5;
Wen = 1055.9;
dw = 0.01;

Copt = struct('Nmax', 300, 'angopt', 1e-1, 'DynDscale', 1);
fmuls = [1 10 20];
U0 = HARMONICSTIFFNESS(GM.M, GM.C, GM.K, Wst, h)\Fl*Famp;

UCs = CONTINUE(@(Uw) GM.HBRESFUN(Uw, Fl*Famp, h, Nt), U0, Wst, Wen, dw, Copt);

% %%
Uoh = kron(eye(Nhc), Rb)*UCs(1:end-1,:);

figure(10)
% clf()
plot(UCs(end,:), vecnorm(Uoh(1:2,:))/Famp, '.-'); hold on
xlim([Wst Wen])