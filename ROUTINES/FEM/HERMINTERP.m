function [Uq, Nint, dNint, qis, eis] = HERMINTERP(X, U, Xq)
%HERMINTERP interpolates the set of values and gradients using Hermite
%Polynomials. Only for 1D functions.
%  
%  USAGE:
%       Xq = HERMINTERP(X, U, Xq);
%  INPUTS:
%       X   : (N,1)
%       U   : (2N,1)
%       Xq  : (Nq,1)
%  OUTPUTS:
%       Uq  : (Nq,1)
    
    U = U(:);
    
    X = X(:);
    Xq = Xq(:);
    
    Les = diff(X);
    Ne  = length(Les);
    Nq = length(Xq);
    
    [qis,eis] = find((X(1:end-1)'-Xq).*(X(2:end)'-Xq)<=0);
    [qis, si] = sort(qis);    eis = eis(si);
    [qis, si] = unique(qis);  eis = eis(si);
    
    if sum(qis(:)'-(1:length(qis)))~=0
%         keyboard
        error('no way')
    end
    
    % Vectorized
    Xes = [X(1:end-1) X(2:end)];
    Les = diff(X);
    
    Xes = Xes(eis, :);
    Les = Les(eis);
    
    xis = (2*Xq-sum(Xes,2))./Les;
    [N, dN] = HERMSF(xis, Les);
    dN = 2./Les.*dN;
    Nint = zeros(Nq, (Ne+1)*2);
    dNint = zeros(Nq, (Ne+1)*2);

    Nint(sub2ind([Nq, (Ne+1)*2], repmat(qis,1,4), (eis-1)*2+(1:4))) = N;
    dNint(sub2ind([Nq, (Ne+1)*2], repmat(qis,1,4), (eis-1)*2+(1:4))) = dN;
    
%     Nint = sparse(Nint);
%     dNint = sparse(dNint);
    
    Uq = Nint*U;
end

function [R, dRdxi] = RESFUN(xi, Xe, Xq, Le)
    [Ns, dNs] = HERMSF(xi, Le);
    
    R = Ns*Xe-Xq;
    dRdxi = dNs*Xe;
end