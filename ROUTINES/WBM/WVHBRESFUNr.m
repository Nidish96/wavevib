function [R, dRdA, dRdw] = WVHBRESFUNr(Ariw, Famp, h, pcs, bcs, joints, Klib)
%WVHBRESFUN returns the residue and Jacobians under Harmonic Balance for
%the Wave-Based Model under reduced parameterization. See "WVAMATr.m" for
%details on reduced parameterization.
%
%   USAGE:
%       [R, dRdA, dRdw] = WVHBRESFUNr(Ariw, h, pcs, bcs, joints, Klib);
%   INPUTS:
%       Ariw    : (Npcs*Nwc*Nhc+1, 1) vector of unknowns + frequency
%       h       : (1, Nh) list of harmonics
%       pcs     : (array of structs) list of wave-based pieces. 
%           See "WVAMAT.m" for documentation.
%       bcs     : (array of structs) list of boundary conditions.
%           See "WVAMAT.m" for documentation.
%       joints  : (array of structs) list of joints.
%           See "WVAMAT.m" for documentation.
%       Klib    : (array of structs) list of dispersion relationships.
%           See "WVAMAT.m" for documentation.
%   OUTPUTS:
%       R       : (Npcs*Nwc*Nhc,1) Residue vector
%       dRdA    : (Npcs*Nwc*Nhc,Npcs*Nwc*Nhc) Jacobian matrix
%       dRdw    : (Npcs*Nwc*Nhc,1) Jacobian wrt frequency

    Nwc = size(pcs(1).wcomps,1);
    Npcs = length(pcs);
    Nhc = sum((h==0)+2*(h~=0));  % Number of Harmonic Terms

    % Linear Part
    w = Ariw(end);
    [Amat, dAmatdw, ~, Fv, dFvdw, ~, JEV] = WVAMATr([w;0], h, pcs, bcs, joints, Klib, 'r');

    % Nonlinearities
    nlis = find(arrayfun(@(j) ~isempty(j.nl), joints));
    nnl = length(nlis);
    FNL = zeros(Npcs*Nwc*Nhc,1);
    dFNLdA = zeros(Npcs*Nwc*Nhc);
    dFNLdw = zeros(Npcs*Nwc*Nhc,1);

    for i=1:nnl
        k = nlis(i);

        U = JEV(i).Lj*Ariw(1:end-1) - JEV(i).LjR*Famp;
        dUdw = JEV(i).dLjdw*Ariw(1:end-1) - JEV(i).dLjRdw*Famp;

        [Fnlk, dFnldUk, dFnldwk] = joints(k).nl([2*U; w]);
        Fnlk = 1/2*Fnlk;
        dFnldUk = 1/2*dFnldUk*2;
        dFnldwk = 1/2*dFnldwk;

        FNL = FNL + JEV(i).Gj*Fnlk;
        dFNLdA = dFNLdA + JEV(i).Gj*dFnldUk*JEV(i).Lj;
        dFNLdw = dFNLdw + JEV(i).dGjdw*Fnlk + ...
            JEV(i).Gj*dFnldUk*dUdw + ...
            JEV(i).Gj*dFnldwk;
    end

    % Setup Residue
    R = Amat*Ariw(1:end-1) + FNL - Fv*Famp;
    dRdA = Amat + dFNLdA;
    dRdw = dAmatdw*Ariw(1:end-1) + dFNLdw - dFvdw*Famp;
end