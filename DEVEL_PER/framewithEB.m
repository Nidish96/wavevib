clc
clear all
addpath('../ROUTINES/SOLVERS/')
addpath('../ROUTINES/WBM/')
%% Setup Model
Ey = 190e9;
rho = 7680;
wid = 0.2;  % Width
brd = 0.4;  % Breadth
Ar = wid*brd;  % Area
Iy = wid^3*brd/12;  % 2nd moment of area
Klib = [struct('K', @(w,xi) sqrt(w)*(rho*Ar/Ey/Iy)^(0.25));
        struct('K', @(w,xi) w*sqrt(Ey/rho))];
%Wave components
wcomps = [1 1;  
         -1 1;  
         1j 1; 
        -1j 1;
         1j 2;
        -1j 2]; 
%pieces
L = 4.0;  
H = 3.0;
pcs = [struct('coords',[0 0;0 2.5;0 H],'wcomps', wcomps);
       struct('coords',[0 H;L H],'wcomps', wcomps);   
       struct('coords',[L H;L 0],'wcomps', wcomps)]; 
% BCs (Fix_Fix)
bcs = [struct('i', 1, 'cofs', @(w,xi) [1 1 1 1 0 0; 1 -1 1j -1j 0 0; 0 0 0 0 1 1]);
       struct('i', 7, 'cofs', @(w,xi) [1 1 1 1 0 0; 1 -1 1j -1j 0 0; 0 0 0 0 1 1])];

%joint
cofs = @(w,xi)[[1 1 1 1 0 0 0 0 0 0 0 0] + (0.5*brd*Klib(1).K(w,xi))*[1 -1 1j -1j 0 0 0 0 0 0 0 0] + [0 0 0 0 0 0 0 0 0 0 1 1];
    [0 0 0 0 0 0 1 1 1 1 0 0] + (0.5*brd*Klib(1).K(w,xi))*[0 0 0 0 0 0 1 -1 1j -1j 0 0] - [0 0 0 0 1 1 0 0 0 0 0 0];
    [1 -1 1j -1j 0 0 -1 1 -1j 1j 0 0];
    (-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0 0 0 0 0 0 0] + (Ey*Ar*Klib(2).K(w,xi))*[0 0 0 0 0 0 0 0 0 0 1j -1j];
    (-Ey*Iy*Klib(1).K(w,xi)^3)*[0 0 0 0 0 0 1 -1 -1j 1j 0 0] - (Ey*Ar*Klib(2).K(w,xi))*[0 0 0 0 1j -1j 0 0 0 0 0 0]; 
    (-Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0 -1 -1 1 1 0 0] + (-0.5*brd*Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0 1 -1 -1j 1j 0 0]];

joints = [struct('type', 2, 'i', 3, 'j', 4, 'cofs', cofs);
          struct('type', 2, 'i', 5, 'j', 6, 'cofs', cofs)];

%% Excitation
Mx = @(w,xi)inv([(-Ey*Iy*Klib(1).K(w,xi)^3)*[1 -1 -1j 1j 0 0];
                 (-Ey*Iy*Klib(1).K(w,xi)^2)*[1 1 -1 -1 0 0];
                  (Ey*Ar*Klib(2).K(w,xi))*[0 0 0 0 1j -1j];
                  [1 1 1 1 0 0];
                  Klib(1).K(w,xi)*[1 -1 1j -1j 0 0];
                  [0 0 0 0 1 1]]);

excs = struct('i', 2, 'nh', 1, 'rcofs', @(w,xi) Mx(w,xi)*[1/2;0;0;0;0;0], ...
    'rcofs0', [1/2;0;0;0;0;0]);

%% Preprocess Everything
[pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib);
Nwc = size(wcomps,1);  % Number of wave components

%% Conduct Linear Forced Response Analysis
Nw = 500;
Ws = linspace(1, 1e3, Nw);
Npts = pcs(end).irange(end);  
ACs = zeros(Npts*Nwc,Nw);
Famp = 200;
for iw=1:Nw
    [Amat, ~, ~, Fv] = WVAMAT([Ws(iw);0], 1, pcs, bcs, joints, Klib);
    ACs(:,iw) = Amat\(Fv*Famp);
end
%% Plot Forced Response
opi = 13:18;  % Output wave coefficients
figure(1)
clf()
subplot(2,1,1)
semilogy(Ws/1e3, abs(2*sum(ACs(opi,:))));
grid on
ylabel('Response (m)')
subplot(2,1,2)
plot(Ws/1e3, rad2deg(angle(2*sum(ACs(opi,:)))));
grid on
set(gca, 'YTick', -180:90:180)
ylabel('Phase (degs)')
xlabel('Frequency (k rad/s)')