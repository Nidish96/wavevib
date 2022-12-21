function [R, dRdU, dRdw, z] = MDOF3D_NLHYST_HBRESFUN(Uw, Pars, L, pA, MESH, M, C, K, Fl, h, Nt, varargin)
%MDOF_NLHYST_HBRESFUN
%
%
  if length(varargin)>=1
    nldofs = varargin{1};
  else 
    nldofs = 1:size(L,1);
  end
  Nhc = sum(h==0) + 2*sum(h~=0);
  Nph = size(L, 1);  % Physical DoFs
  Nd = size(L, 2);   % Projected DoFs

  Ncdof = MESH.Nn*MESH.dpn;  % Number of Contact DoFs
  Nqp = MESH.Ne*MESH.Nq^2;
  Uhqp = zeros(Nqp*3, Nhc, 'single');
  Uhqp(1:3:end, :) = MESH.INTERP_QP(L(1:3:Ncdof, :)*reshape(Uw(1:end-1), Nd, Nhc));  % Qm*L*Ux; (Nqp x Nhc)
  Uhqp(2:3:end, :) = MESH.INTERP_QP(L(2:3:Ncdof, :)*reshape(Uw(1:end-1), Nd, Nhc));  % Qm*L*Uy; (Nqp x Nhc)
  Uhqp(3:3:end, :) = MESH.INTERP_QP(L(3:3:Ncdof, :)*reshape(Uw(1:end-1), Nd, Nhc));  % Qm*L*Uz; (Nqp x Nh)

  % NONLINEARITY
  [FQP, JQP] = MESH.fcont(reshape(Uhqp, 3*Nqp*Nhc, 1), h, Nt, pA*Pars);
  FQP = double(FQP);
  JQP = double(JQP);
  
  % % Quadrature Integration
  % Fnl = zeros(Nph*Nhc,1, 'single');
  % Jnl = zeros(Nph*Nhc, 'single');
  % for ih=1:Nhc
  %   for fi=1:3    
  %     Fnl((ih-1)*Nph+(fi:3:Ncdof)) = MESH.Tm*FQP((ih-1)*3*Nqp+(fi:3:3*Nqp));  % need to typecast for operation to work

  %     for jh=1:Nhc
  %       for fj=1:3
  %         Jnl((ih-1)*Nph+(fi:3:Ncdof), (jh-1)*Nph+(fj:3:Ncdof)) = ...
  %         MESH.Tm*JQP((ih-1)*3*Nqp+(fi:3:3*Nqp), (ih-1)*3*Nqp+(fj:3:3*Nqp))*MESH.Qm;
  %       end
  %     end
  %   end
  % end
  % % Apply Projection "L"
  % if Nph>Nd  % L represents some subspace restriction, so the following will work
  %   for ih=1:Nhc
  %     Fnl((ih-1)*Nd+(1:Nd)) = L'*Fnl((ih-1)*Nph+(1:Nph));
      
  %     Jnl((ih-1)*Nd+(1:Nd), :) = L'*Jnl((ih-1)*Nph+(1:Nph), :);
  %     Jnl(:, (ih-1)*Nd+(1:Nd)) = Jnl(:, (ih-1)*Nph+(1:Nph))*L;
  %   end

  %   Fnl = Fnl(1:Nd*Nhc);
  %   Jnl = Jnl(1:Nd*Nhc, 1:Nd*Nhc);
  % end
  
  
  Fnl = zeros(Nd*Nhc,1, 'single');
  Jnl = zeros(Nd*Nhc, 'single');  
  for ih=1:Nhc
      for fi=1:3
          Fnl((ih-1)*Nd+(1:Nd)) = Fnl((ih-1)*Nd+(1:Nd)) + ...
              L(fi:3:Ncdof,:)'*MESH.Tm*FQP((ih-1)*3*Nqp+(fi:3:3*Nqp));
          
          for jh=1:Nhc
              for fj=1:3
                  Jnl((ih-1)*Nd+(1:Nd), (jh-1)*Nd+(1:Nd)) = Jnl((ih-1)*Nd+(1:Nd), (jh-1)*Nd+(1:Nd)) + ... 
                      L(fi:3:Ncdof,:)'*MESH.Tm*JQP((ih-1)*3*Nqp+ ...
                                                   (fi:3:3*Nqp), ...
                                                   (jh-1)*3*Nqp+(fj:3:3*Nqp))*MESH.Qm*L(fj:3:Ncdof, :);
              end
          end
      end
  end

				% Residue
  [E, dEdw] = HARMONICSTIFFNESS(M, C, K, Uw(end), h);
  
  R = E*double(Uw(1:end-1)) + Fnl - Fl;
  dRdU = single(full(E)) + Jnl;
  dRdw = single(dEdw*double(Uw(1:end-1)));
end
