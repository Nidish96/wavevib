function [zinds,hinds,rinds0,rinds,iinds] = HINDS(Ndofs, h)
%HINDS produces the indices required for converting complex representation
%to real-imaginary representation in the presence of multiple DOFs.
%
%   USAGE:
%       [rinds0,zinds,hinds,rinds,iinds] = HINDS(Ndofs, h);
%   INPUTS:
%       Ndofs   : (int) 
%       h       : (1,Nh)
%   OUTPUTS:
%       zinds   : indices of zero harmonics in complex rep
%       hinds   : indices of non-zero harmonics in complex rep
%       rinds0  : indices where zero harmonics will go in RI rep
%       rinds   : indices where real parts will go in RI rep
%       iinds   : indices where imag parts will go in RI rep
    Nh = length(h);
    if h(1)==0
        rinds0 = 1:Ndofs;
        zinds = 1:Ndofs;
        hinds = (Ndofs+1):(Ndofs*Nh);
        rinds = repmat(1:Ndofs, 1, Nh-1) + kron((0:Nh-2)*2*Ndofs+Ndofs, ones(1,Ndofs));  % real indices
        iinds = repmat(1:Ndofs, 1, Nh-1) + kron((1:Nh-1)*2*Ndofs, ones(1,Ndofs));  % imag indices
    else
        rinds0 = [];
        zinds = [];
        hinds = 1:Ndofs*Nh;
        rinds = repmat(1:Ndofs, 1, Nh) + kron((0:Nh-1)*2*Ndofs, ones(1,Ndofs));  % real indices
        iinds = repmat(1:Ndofs, 1, Nh) + kron((0:Nh-1)*2*Ndofs+Ndofs, ones(1,Ndofs));  % imag indices
    end
end