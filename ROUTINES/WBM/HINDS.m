function [zinds,hinds,rinds0,rinds,iinds] = HINDS(Ndofs, h)
%HINDS produces the indices required for converting complex representation
%to real-imaginary representation in the presence of multiple DOFs.
%   Complex representation: [A1; A2; A3; ...]
%     where each "Ai" contains the values of all the DOFs for harmonic "i".
%   Real-Imaginary representation: [real(A1);imag(A1);real(A2);imag(A2);..]
%
%       The persence of zero harmonics is handed specially (no imaginary
%               coefficients).
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
%         rinds = repmat(1:Ndofs, 1, Nh-1) + kron((0:Nh-2)*2*Ndofs+Ndofs, ones(1,Ndofs));  % real indices
%         iinds = repmat(1:Ndofs, 1, Nh-1) + kron((1:Nh-1)*2*Ndofs, ones(1,Ndofs));  % imag indices

        rinds = reshape((1:Ndofs)' + (0:Nh-2)*2*Ndofs+Ndofs, 1, []);
        iinds = rinds + Ndofs;
    else
        rinds0 = [];
        zinds = [];
        hinds = 1:Ndofs*Nh;
%         rinds = repmat(1:Ndofs, 1, Nh) + kron((0:Nh-1)*2*Ndofs, ones(1,Ndofs));  % real indices
%         iinds = repmat(1:Ndofs, 1, Nh) + kron((0:Nh-1)*2*Ndofs+Ndofs, ones(1,Ndofs));  % imag indices

        rinds = reshape((1:Ndofs)' + (0:Nh-1)*2*Ndofs, 1, []);
        iinds = rinds + Ndofs;
    end
end