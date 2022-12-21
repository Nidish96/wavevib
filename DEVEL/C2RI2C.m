function [Uout, L] = C2RI2C(Uin, h, drc)
%C2RI2C converts representation of list of harmonics in complex to
%real-imaginary form and vise versa. The zero harmonic term (if present) is
%handled appropriately.
%
%   USAGE:
%       [Uout] = C2RI2C(Uin, h, dir);
%   INPUTS:
%       Uin     : (Ndof*Nh, Npts) 
%                   OR
%                 (Ndof*Nhc, Npts)  
%       h       : (1, Nh)
%       drc     : 'c2ri' or 'ri2c'
%   OUTPUTS:
%       Uout    : (Ndof*Nhc, Npts)
%                   OR
%                 (Ndof*Nh, Npts)
%       L       : (Ndof*Nh, Ndof*Nhc)  % ri2c Tfmx

    
    Nh = length(h);
    Nhc = sum((h==0)+2*(h~=0));

    N = size(Uin,1);
    Npts = size(Uin,2);
    drc = lower(drc);
    if strcmp(drc, 'c2ri')
        if mod(N, Nh)~=0
            error('Unrecognized number of terms provided.');
        end
        Ndofs = N/Nh;
        rinds = repmat(1:Ndofs, 1, Nh) + kron(((1:Nh)-1)*2*Ndofs, ones(1,Ndofs));  % real indices
        iinds = repmat(1:Ndofs, 1, Nh) + kron(((1:Nh)-1)*2*Ndofs+Ndofs, ones(1,Ndofs));  % imag indices

        Uout = zeros(2*N, Npts);
        Uout(rinds,:) = real(Uin);
        Uout(iinds,:) = imag(Uin);
        if h(1)==0
            Uout(Ndofs+(1:Ndofs),:) = [];
        end
    elseif strcmp(drc, 'ri2c')
        if mod(N, Nhc)~=0
            error('Unrecognized number of terms provided.');
        end
        Ndofs = N/Nhc;
        rinds = repmat(1:Ndofs, 1, Nh) + kron(((1:Nh)-1)*2*Ndofs, ones(1,Ndofs));  % real indices
        iinds = repmat(1:Ndofs, 1, Nh) + kron(((1:Nh)-1)*2*Ndofs+Ndofs, ones(1,Ndofs));  % imag indices
        
        if h(1)==0
            Uin = [Uin(1:Ndofs,:); zeros(Ndofs,Npts);Uin(Ndofs+1:end,:)];
        end
        Uout = Uin(rinds,:)+1j*Uin(iinds,:);
    else
        error('Unrecognized drc option.');
    end

    % Transformation Matrix
    rinds = repmat(1:Ndofs, 1, Nh) + kron(((1:Nh)-1)*2*Ndofs, ones(1,Ndofs));  % real indices
    iinds = repmat(1:Ndofs, 1, Nh) + kron(((1:Nh)-1)*2*Ndofs+Ndofs, ones(1,Ndofs));  % imag indices

    L = zeros(Nh*Ndofs, 2*Nh*Ndofs);
    L(:, rinds) = 1;
    L(:, iinds) = 1j;
    if h(1)==0
        L(:, Ndofs+(1:Ndofs)) = [];
    end
end