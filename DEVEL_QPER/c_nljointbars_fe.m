clc
clear all

%% Parameters
Ey = 2.62e11;
rho = 1280;
Ar = 31.75e-3*54e-3;

wcomps = [1j 1; -1j 1];

ell = 1;
xF1 = 0.28;
xF2 = 4/6;

% Joint Properties
xJ = 1/3;

kJ = 1e9;
cJ = 320;
gJ = 1e8;

%% Finite Element Matrices
Me = @(Le) rho*Ar*Le/6*[2 1;1 2];
Ke = @(Le) Ey*Ar/Le*[1 -1;-1 1];

Ne1 = 10;
Xn1 = linspace(0, xJ, Ne1+1);
[~, iF1] = min(abs(Xn1-xF1));  iF1 = iF1(1);  % Forcing node
Xn1(iF1) = xF1;
M1 = zeros(Ne1+1);
K1 = zeros(Ne1+1);
for ei=1:Ne1
    M1(ei:ei+1, ei:ei+1) = M1(ei:ei+1, ei:ei+1) + Me(diff(Xn1(ei:ei+1)));
    K1(ei:ei+1, ei:ei+1) = K1(ei:ei+1, ei:ei+1) + Ke(diff(Xn1(ei:ei+1)));    
end

Ne2 = Ne1*2;
Xn2 = linspace(xJ, ell, Ne2+1);
[~, iF2] = min(abs(Xn1-xF2));  iF2 = iF2(1);  % Forcing node
Xn1(iF2) = xF2;
M2 = zeros(Ne2+1);
K2 = zeros(Ne2+1);
for ei=1:Ne2
    M2(ei:ei+1, ei:ei+1) = M2(ei:ei+1, ei:ei+1) + Me(diff(Xn2(ei:ei+1)));
    K2(ei:ei+1, ei:ei+1) = K2(ei:ei+1, ei:ei+1) + Ke(diff(Xn2(ei:ei+1)));    
end
iF2 = Ne1+1+iF2;

M = blkdiag(M1, M2);
K = blkdiag(K1, K2);
Lb = eye(Ne1+Ne2+2);
Lb = Lb(:, 2:end-1);
Lnl = Lb(Ne1+2,:)-Lb(Ne1+1,:);
Fvb1 = Lb(iF1, :)';
Fvb2 = Lb(iF2, :)';
Rb = Lb(Ne1+2,:);
Mb = Lb'*M*Lb;
Kb = Lb'*K*Lb;

KJ = Lnl'*kJ*Lnl;
CJ = Lnl'*cJ*Lnl;

%% Setup ode45 solve
ws = 78e3*[1;pi];
F0 = 150e5;

fex = @(t, i) F0*cos(ws(i)*t);
fun = @(t, y) [zeros(Ne1+Ne2) eye(Ne1+Ne2); -Mb\(Kb+KJ) -Mb\CJ]*y - [zeros(Ne1+Ne2,1); Lnl']*gJ*(Lnl*y(1:Ne1+Ne2))^3 + [zeros(Ne1+Ne2,1); Mb\(Fvb1*fex(t,1)+Fvb2*fex(t,2))];

T0 = 0;
T1 = 1e-1;
y0 = zeros((Ne1+Ne2)*2,1);

[T, Y] = ode45(fun, [T0 T1], y0);
y = Rb*Y(:, 1:Ne1+Ne2)';
save('TransFEres.mat', 'T', 'y')
%% Plot Results
figure(1)
clf()
plot(T, y)

xlabel('Time (s)')
ylabel('Displacement (m)')