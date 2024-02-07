function [R, dRdA, dRdw] = WVHBRESFUN(Ariw, Famp, h, pcs, bcs, joints, Klib)
%WVHBRESFUN returns the residue and Jacobians under Harmonic Balance for
%the Wave-Based Model.
%       This works for the Quasi-Periodic case.
%
%   USAGE:
%       [R, dRdA, dRdw] = WVHBRESFUN(Ariw, h, pcs, bcs, joints, Klib);
%   INPUTS:
%       Ariw    : (Npts*Nwc*Nhc+Nc, 1) vector of unknowns + frequency
%       h       : (Nh, Nc) list of harmonics
%       pcs     : (array of structs) list of wave-based pieces. 
%           See "WVAMAT.m" for documentation.
%       bcs     : (array of structs) list of boundary conditions.
%           See "WVAMAT.m" for documentation.
%       joints  : (array of structs) list of joints.
%           See "WVAMAT.m" for documentation.
%       Klib    : (array of structs) list of dispersion relationships.
%           See "WVAMAT.m" for documentation.
%   OUTPUTS:
%       R       : (Npts*Nwc*Nhc,1) Residue vector
%       dRdA    : (Npts*Nwc*Nhc,Npts*Nwc*Nhc) Jacobian matrix
%       dRdw    : (Npts*Nwc*Nhc,Nc) Jacobian wrt frequency

    Nwc = size(pcs(1).wcomps,1);
    assert(all(arrayfun(@(p) size(p.wcomps,1), pcs(2:end))==Nwc));
    Npts = pcs(end).irange(end);
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
    Nc = size(h,2);
    
    % Linear Parts
    ws = Ariw(end-Nc+1:end);
    [Amat, dAmatdw, ~, Fv, dFvdw, ~, JEV] = WVAMAT([ws;0], h, pcs, bcs, joints, Klib, 'r');

    % Nonlinearities
    nlis = find(arrayfun(@(j) ~isempty(j.nl), joints));
    nnl = length(nlis);
    FNL = zeros(Npts*Nwc*Nhc,1);
    dFNLdA = zeros(Npts*Nwc*Nhc);
    dFNLdw = zeros(Npts*Nwc*Nhc,Nc);

    for i=1:nnl
        k = nlis(i);

        U = JEV(i).Lj*Ariw(1:end-Nc);
        dUdw = cell2mat(arrayfun(@(a) JEV(i).dLjdw(:, :, a)*Ariw(1:end-Nc), 1:Nc, 'UniformOutput', false));

        nfnl = size(JEV(i).Gj,2)/Nhc;
        Fsc = 1/2*ones(nfnl*Nhc,1);
        Usc = 2*ones(joints(k).nld*Nhc,1);
        if sum(all(h==0,2))~=0
            assert(all(h(1,:)==0))
            Fsc(1:nfnl) = 2*Fsc(1:nfnl);
            Usc(1:joints(k).nld) = Usc(1:joints(k).nld)/2;
        end
        
        [Fnlk, dFnldUk, dFnldwk] = joints(k).nl([Usc.*U; ws]);
        
        Fnlk = Fsc.*Fnlk;
        dFnldUk = Fsc.*dFnldUk.*Usc';
        dFnldwk = Fsc.*dFnldwk;

        FNL = FNL + JEV(i).Gj*Fnlk;
        dFNLdA = dFNLdA + JEV(i).Gj*dFnldUk*JEV(i).Lj;
        dFNLdw = dFNLdw + cell2mat(arrayfun(@(a) JEV(i).dGjdw(:, :, a)*Fnlk, 1:Nc, 'UniformOutput', false)) + ...
            JEV(i).Gj*dFnldUk*dUdw + ...
            JEV(i).Gj*dFnldwk;
    end

    % Setup Residue
    R = Amat*Ariw(1:end-Nc) + FNL - Fv*Famp;
    dRdA = Amat + dFNLdA;
    dRdw = cell2mat(arrayfun(@(a) dAmatdw(:, :, a)*Ariw(1:end-Nc), 1:Nc, 'UniformOutput', false)) + dFNLdw - dFvdw*Famp;
end
