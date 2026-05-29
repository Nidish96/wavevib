function [FNL, dFNLdU, dFNLdw] = HFRIC(Uw, mu, kt, kn, gap, h, Nt)
%HFRIC returns the Fourier Coefficients of the fricitonal contact model 
% implemented as a linear-penalty model in the normal direction and  
% a jenkin's model in the tangential direction. 
% Force that is returned is:
%                   fn = kn*max(un-gap, 0)
%                   ft = kt*(ut-ut0) + ft0, stick
%                        mu*fn*signm(fstick), slip
%                        0, separation
%   positive (u-gap) quantity is considered as contact and negative is
%   considered as separation (zero force). 
%
%   USAGE: 
%       [FNL, dFNLdU, dFNLdw] = HFRIC(Uw, mu, kt, kn, gap, h, Nt);
%   INPUTS:
%       Uw            : (Nhc*2+1,1) Vector of nonlinear relative DOFs+freq.
%                       Convention -> [a0t a0n a1t a1n b1t b1n a2t ..., freq]^T 
%       mu,kt,kn,gap  : (scalar) or (Nd,1) parameters]
%       h             : (Nh,1) Vector of Harmonics
%       Nt            : (int) Number of samples for AFT
%   OUTPUTS:
%       FNL     : (Nhc*2, 1) Force harmonics
%       dFNLdU  : (Nhc*2,Nhc*2) Force harmonic jacobian wrt U
%       dFNLdw  : (Nhc*2,1) Force harmonic jacobian wrt w
    
    Nhc = sum((h==0)+2*(h~=0));
    cst = AFT(eye(Nhc), h, Nt, 'f2t');

    Nd = (size(Uw,1)-1)/Nhc;  % Number of DOFs
    assert(Nd==2);

    % Time domain displacements and force declaration
    utn = AFT(reshape(Uw(1:end-1), 2,Nhc)', h,Nt, 'f2t');  % (Nt, 2)
    ftn = zeros(Nt, 2);

    % Normal force
    ftn(:, 2) = kn*max(utn(:,2)-gap, 0);
    dfndun = kn*(utn(:,2)>gap);
    dfndUn = dfndun.*cst;

    % Tangential force
    dftdUn = zeros(Nt, Nhc);
    dftdUt = zeros(Nt, Nhc);
    for itn = 1:2
        for ti = 1:Nt
            tim1 = mod(ti-1 -1,Nt)+1;
    
            fsp = kt*(utn(ti,1)-utn(tim1,1)) + ftn(tim1,1); % stick prediction
            if abs(fsp)<mu*ftn(ti,2)  % Stuck
                ftn(ti,1) = fsp;
                dftdUt(ti, :) = kt*(cst(ti,:)-cst(tim1,:)) + dftdUt(tim1, :);
                dftdUn(ti, :) = dftdUn(tim1, :);
            else  % Slipped
                ftn(ti,1) = mu*ftn(ti,2)*sign(fsp);
                dftdUt(ti, :) = 0;
                dftdUn(ti, :) = mu*dfndUn(ti, :)*sign(fsp);
            end
        end
    end

    tmp = AFT([ftn dftdUt dftdUn dfndUn], h,Nt, 't2f');
    Ftn = tmp(:, 1:2);
    dFtdUt = tmp(:, 2+(1:Nhc));
    dFtdUn = tmp(:, 2+Nhc+(1:Nhc));
    dFndUn = tmp(:, 2+2*Nhc+1:end);

    % AFT for Force 
    FNL = reshape(Ftn', 2*Nhc, 1);
    dFNLdU = zeros(2*Nhc, 2*Nhc);
    dFNLdU(1:2:end, 1:2:end) = dFtdUt;
    dFNLdU(1:2:end, 2:2:end) = dFtdUn;
    dFNLdU(2:2:end, 2:2:end) = dFndUn;
    dFNLdw = zeros(Nhc*2,1);
end