function [P_nat] = TRI2D_MAP2NAT(P_glob,V)
%TRIE2D_MAP2NAT Returns the input points (in global CS) in
%equivalent natural CS
% USAGE:
%	[P_nat] = TRI2D_MAP2NAT(P_glob,V);
% INPUTS:
%   P_glob	: Npx2 query points (global CS)
%   V 		: 4x2 ordered coordinates of vertices (global CS)
% OUTPUTS:
%   P_nat	: Npx2 mapped points (natural CS)

    opt = optimoptions('fsolve','specifyObjectiveGradient',true, ...
                       'display','off');
    Np = size(P_glob,1);
    P_nat = zeros(Np,2);
%     for n=1:Np
%         P_nat(n,:) = fsolve(@(x) RESFUN(x,V,P_glob(n,:)), P_nat(n,:), ...
%                             opt); 
%     end
    P_nat = fsolve(@(x) RESFUN(x,V,P_glob),P_nat,opt);
    P_nat(P_nat>1) = 1.0;
    P_nat(P_nat<0) = 0.0;
end

function [R,J] = RESFUN(x,V,P_glob)
%RESFUN Residual function for mapping
% USAGE:
%	[R,J] = RESFUN(x,V,P_glob);
% INPUTS:
%   x,V,P_glob 
% OUTPUTS:
%   R,J

    Np = size(P_glob,1);
    R = TRI2D_SF(x)*V-P_glob;
    J = TRI2D_SD(x)*V;
    J = (kron(eye(Np),ones(2))).*kron(ones(1,Np),J);
end