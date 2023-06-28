function [P_nat] = QUAD2D_MAP2NAT(P_glob,V)
%QUAD2D_MAP2NAT Returns the input points (in global CS) in
%equivalent natural CS
% USAGE:
%	[P_nat] = QUAD2D_MAP2NAT(P_glob,V);
% INPUTS:
%   P_glob	: Npx2 query points (global CS)
%   V 		: 4x2 ordered coordinates of vertices (global CS)
% OUTPUTS:
%   P_nat	: Npx2 mapped points (natural CS)

    opt = optimoptions('fsolve','specifyObjectiveGradient',false, ...
                       'display','off','tolfun',1e-12,'FunctionTolerance',1e-15);
%     opt = optimoptions('fmincon','display','off');
    Np = size(P_glob,1);
    P_nat = zeros(Np,2);
    for n=1:Np
        [P_nat(n,:),~,eflg] = fsolve(@(x) RESFUN(x,V,P_glob(n,:)), P_nat(n,:), ...
                                opt);
        if eflg<=0
            func = @(x) vecnorm(QUAD2D_SF(x)*V-P_glob)/vecnorm(P_glob);
            P_nat(n,:) = fmincon(func,[0,0],[],[],...
                [],[],[-1 -1],[1 1],[],optimoptions('fmincon','display','off'));
        end
    end
%     P_nat = fsolve(@(x) RESFUN(x,V,P_glob),P_nat,opt);
%     P_nat(abs(P_nat)>1) = sign(P_nat(abs(P_nat)>1));
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
    R = QUAD2D_SF(x)*V-P_glob;
    J = QUAD2D_SD(x)*V;
%     J = (kron(eye(Np),ones(2))).*kron(ones(1,Np),J);   
end