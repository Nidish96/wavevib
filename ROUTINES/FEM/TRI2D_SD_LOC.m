function [dN] = TRI2D_SD_LOC(P,V)
%TRI2D_SLOC Returns the shape function derivatives of a 2D Tri element 
%evaluated at query local points
% USAGE:
%	[N] = TRI2D_SD_LOC(P,V);
% INPUTS:
%   P 		: Npx2 query points (local CS)
%   V       : 3x2 nodal coordinates
% OUTPUTS:
%   N		: Npx3 shape function evaluates
    
    Np = size(P,1);
    Amat = [ones(3,1) V(:,1) V(:,2)];
    dN = zeros(2*Np,3);
    dN(1:2:end,:) = kron(ones(Np,1),[0 1 0])/Amat;
    dN(2:2:end,:) = kron(ones(Np,1),[0 1 0])/Amat;
end