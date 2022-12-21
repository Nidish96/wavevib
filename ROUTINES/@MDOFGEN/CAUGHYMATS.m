function [MATS] = CAUGHYMATS(m, num, typ)
%CAUGHYMATS returns matrices from the Caughy Series
%
%   USAGE:
%       MATS = m.CAUGHYMATS(num, typ);
%   INPUTS:
%       m       : MDOFGEN Object
%       num     : number of Caughy Matrices
%       typ     : [0,1] 
%                   0: M*(inv(M)*K)^n (M, K, K*inv(M)*K, ...)
%                   1: K*(inv(K)*M)^n (K, M, M*inv(K)*M, ...)
%   OUTPUTS:
%       MATS    : (m.Ndofs, m.Ndofs, num) pages of matrices
    
    MATS = zeros(m.Ndofs, m.Ndofs, num);
    MATS(:, :, 1) = m.M;
    if num>=2
        MATS(:, :, 2) = m.K;
    end
    switch typ
        case 0
            tmp = m.M\m.K;
        case 1
            tmp = m.K\m.M;
        otherwise
            error('Unknown Caughy Series type')
    end
    for ni=3:num
        MATS(:, :, ni) = MATS(:, :, ni-1)*tmp;
    end
end
