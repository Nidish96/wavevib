function [R, dRdA, dRdw] = WVHBRESFUN(Ariw, Famp, h, pcs, bcs, joints, Klib, wcomps)
%WVHBRESFUN
%
%   USAGE:
%       [R, dRdA, dRdw] = WVHBRESFUN(Ariw, h, pcs, bcs, joints, Klib, wcomps);
%
%   INPUTS:
%       
%   OUTPUTS:
    Nwc = size(wcomps,1);
    Npts = pcs(end).irange(end);
    Nhc = sum((h==0)+2*(h~=0));  % Number of Harmonic Terms
    Nh = length(h);

    [zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);
    ac = zeros(Npts*Nwc*Nhc,1);
    ac([rinds0 rinds iinds]) = [real(Ariw([zinds hinds])); imag(Ariw(hinds))];

    % Linear Part
    w = Ariw(end);
    [Amat, dAmatdw, ~, Fv, dFvdw, ~, JEV] = WVAMAT([w;0], h, pcs, bcs, joints, Klib, wcomps, 'r');

    % Nonlinearities
    nlis = find(arrayfun(@(j) ~isempty(j.nl), joints));
    nnl = length(nlis);
    FNL = zeros(Npts*Nwc*Nhc,1);
    dFNLdA = zeros(Npts*Nwc*Nhc);
    dFNLdw = zeros(Npts*Nwc*Nhc,1);
    for i=1:nnl
        k = nlis(i);

        U = JEV.Lj(i:nnl:end,:)*Ariw(1:end-1);
        dUdw = JEV.dLjdw(i:nnl:end,:)*Ariw(1:end-1);

        [Fnlk, dFnldUk, dFnldwk] = joints(k).nl([2*U; w]);
        Fnlk = 1/2*Fnlk;
        dFnldUk = 1/2*dFnldUk*2;
        dFnldwk = 1/2*dFnldwk;

        FNL = FNL + JEV.Gj(:,i:nnl:end)*Fnlk;
        dFNLdA = dFNLdA + JEV.Gj(:,i:nnl:end)*dFnldUk*JEV.Lj(i:nnl:end,:);
        dFNLdw = dFNLdw + JEV.dGjdw(:,i:nnl:end)*Fnlk + ...
            JEV.Gj(:,i:nnl:end)*dFnldUk*dUdw + ...
            JEV.Gj(:,i:nnl:end)*dFnldwk;
    end

    % Setup Residue
    R = Amat*Ariw(1:end-1) + FNL - Fv*Famp;
    dRdA = Amat + dFNLdA;
    dRdw = dAmatdw*Ariw(1:end-1) + dFNLdw - dFvdw*Famp;
end