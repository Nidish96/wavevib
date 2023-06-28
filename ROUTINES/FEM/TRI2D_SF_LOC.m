function [N] = TRI2D_SF_LOC(P,V)
%TRI2D_SF_LOC Returns the shape functions of a 2D Tri element evaluated
%at query local points
% USAGE:
%	[N] = TRI2D_SF_LOC(P,V);
% INPUTS:
%   P 		: Npx2 query points (local CS)
%   V       : 3x2 nodal coordinates
% OUTPUTS:
%   N		: Npx3 shape function evaluates
    
    Amat = [ones(3,1) V(:,1) V(:,2)];
    N = [ones(size(P,1),1) P(:,1) P(:,2)]/Amat;
end