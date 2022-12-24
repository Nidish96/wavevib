function [Det, dDetdw, dDetdxi] = WVLDETFUN(wxi, h, pcs, bcs, joints, Klib)
%WVLDETFUN returns the determinant of the "A" matrix of a wave-based model
%along with its derivatives. Internally, the reduced parameterization is
%employed for computational savings.
%   The real-imaginary representation is employed, so the determinant is
%   always real.
%
%   USAGE:
%       [Det, dDetdw, dDetdxi] = WVLDETFUN(wxi, h, pcs, bcs, joints, Klib);
%   INPUTS:
%       wxi     : (2,1) vector of (frequency,parameter)
%       h       : (1,Nh) list of harmonics
%       pcs     : (array of structs) list of wave-based pieces. 
%           See "WVAMAT.m" for documentation.
%       bcs     : (array of structs) list of boundary conditions.
%           See "WVAMAT.m" for documentation.
%       joints  : (array of structs) list of joints.
%           See "WVAMAT.m" for documentation.
%       Klib    : (array of structs) list of dispersion relationships.
%           See "WVAMAT.m" for documentation.
%   OUTPUTS:
%       Det     : (scalar,real) determinant
%       dDetdw  : (scalar,real) 
%       dDetdxi : (scalar,real)
    [Amat, dAmatdw, dAmatdxi] = WVAMATr(wxi, h, pcs, bcs, joints, Klib, 'r');

    Det = det(Amat);
    tmp = Amat\[dAmatdw dAmatdxi];
    dDetdw = Det*trace(tmp(1:size(Amat,1),1:size(Amat,1)));
    dDetdxi = Det*trace(tmp(1:size(Amat,1),size(Amat,1)+(1:size(Amat,1))));
end