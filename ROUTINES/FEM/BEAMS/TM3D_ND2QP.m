function [Q, T] = TM3D_ND2QP(Les, Wys, Zs, No)
%TM3D_ND2QP produces node to quadrature point interpolation and
%integration matrices.
%
%  USAGE:
%    [Q, T] = TM3D_ND2QP(Les, Wys, No);
%  INPUTS:
%    Les	:
%    Wys	:
%    Zs		: 
%    No		:
%  OUTPUTS:
%    Q		: Interpolation
%    T 		: Integration
  
  if length(Les)~=length(Wys) || length(Les)~=length(Zs)
    error('wth bro')
  end
  
  [xi, wi] = LGWT(No, -1, 1);
  [xi, si] = sort(xi);
  wi = wi(si);

  Ne = length(Les);
  Q  = zeros(Ne*No^2*3, (Ne+1)*6);
  T  = zeros((Ne+1)*6, Ne*No^2*3);
  for e = 1:Ne
    NS = SFUNCS(xi);

    Ys = NS*[-0.5; 0.5]*Wys(e);
    interp_mx = zeros(No^2*3, 12);
    integ_mx  = zeros(12, No^2*3);
    for yi=1:No
      interp_mx((yi-1)*No*3+(1:No*3), :) = kron(NS, [eye(3), [0 Zs(e) -Ys(yi);-Zs(e) 0 0;Ys(yi) 0 0]]);
%       interp_mx((yi-1)*No*3+(1:No*3), :) = kron(NS, [eye(3), -[0 Zs(e) -Ys(yi);-Zs(e) 0 0;Ys(yi) 0 0]]);      

      integ_mx(:, (yi-1)*No*3+(1:No*3)) = kron(Les(e)/2*wi.*NS, [eye(3), [0 Zs(e) -Ys(yi);-Zs(e) 0 0;Ys(yi) 0 0]])'*wi(yi)*Wys(e)/2;
%       integ_mx(:, (yi-1)*No*3+(1:No*3)) = kron(Les(e)/2*wi.*NS, [eye(3), -[0 Zs(e) -Ys(yi);-Zs(e) 0 0;Ys(yi) 0 0]])'*wi(yi)*Wys(e)/2;
    end

    Q((e-1)*(No^2*3)+(1:No^2*3), (e-1)*6+(1:12)) = ...
    Q((e-1)*(No^2*3)+(1:No^2*3), (e-1)*6+(1:12)) + interp_mx;

    T((e-1)*6+(1:12), (e-1)*(No^2*3)+(1:No^2*3)) = ...
    T((e-1)*6+(1:12), (e-1)*(No^2*3)+(1:No^2*3)) + integ_mx;
    
  end  
end

function [Ns] = SFUNCS(xi)
  xi = xi(:);
  Np = length(xi);
  Ns = zeros(Np, 2);
  Ns = [(1-xi)/2, (1+xi)/2];
end
