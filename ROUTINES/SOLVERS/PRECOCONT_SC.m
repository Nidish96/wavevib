function [Ulc, dUlds, Ss] = PRECOCONT_SC(func, u0, lstart, lend, ds, varargin)
%%PRECONT_SC Applies the solver directly for continuation
%
%   USAGE:
%       [Uc] = PRECOCONT_SC(func, u0, l0, l1, dl, varargin)
%   INPUTS:
%       func, u0, l0, l1, dl, Copt
%   OUTPUTS:
%       Uc

  %% Default options
  Copt = struct('Nmax', 100, 'dsmax', ds*5, 'dsmin', ds/5, 'stepadapt', 1, ...
                'angopt', pi/6, 'startdir', sign(lend-lstart),...
                'Display', true, 'itopt', 10, ...
                'itDisplay', false, 'lsrch', 0, 'lsit', 10, ...
                'predictororder', 1, ...  % 0, 1, 3, 5, ...
                'arclengthparm', 'orthogonal',...  % 'orthogonal', 'arclength'
                'Dscale', ones(length(u0)+1,1), ...  % Parameter weight
                'ITMAX', 20, 'DynScale', 0, ...
                'opts', struct('reletol', 1e-12, 'rtol', 1e-12, 'etol', ...
                              1e-12, 'utol', 1e-12, 'crit', 14, ...
                              'Dscale', ones(length(u0),1), ...
                              'Display', false));
  if nargin==6
    nflds = fieldnames(varargin{1});
    for i=1:length(nflds)
      Copt.(nflds{i}) = varargin{1}.(nflds{i});
    end
  end
  Copt.opts.lsrch = Copt.lsrch;
  Copt.opts.lsit = Copt.lsit;
  Copt.opts.Display = Copt.itDisplay;
  Copt.opts.ITMAX = Copt.ITMAX;
  
  %% 
  if Copt.predictororder>3
      warning('High order predictors are known to perform badly around turning points')
  end
  
  %% Allocations
  Ss = zeros(Copt.Nmax, 1);
  Ulc = zeros(length(u0)+1, Copt.Nmax);
  dUlds = zeros(length(u0)+1, Copt.Nmax);

  ds0 = ds;
  
  %% Initial Point Correction
  Copt.opts.Display = true;
  [Ulc(1:end-1,1), ~, eflag] = GSOLVE(@(u) func([u; lstart]), u0, Copt.opts);
  Copt.opts.Display = Copt.itDisplay;
  Ulc(end, 1) = lstart;
  if eflag < 0
      error('Initial point non-convergent!');
  elseif Copt.Display
      disp('Initial Point Converged')
      fprintf('Starting Continuation from %f to %f\n', lstart, lend);
  end
  
  %% Initial Tangent (Unscaled Coordinates)
  [R, dRdU, dRdl] = func(Ulc(:,1));
  dRdU = dRdU*diag(Copt.Dscale(1:end-1)); dRdl = dRdl*Copt.Dscale(end);  
  z  = -dRdU\dRdl;
  al = Copt.startdir/sqrt(1+z'*z);
  dUlds(:, 1) = [z; 1]*al;
  
  %% BEGIN CONTINUATION
  l = lstart;
  n = 1;
  
  ds_ref = ds;
  if Copt.Display
      fprintf('N Lam ds Its: Flag\n');
  end
  Copt.opts.lsrch = 0;  % Line Search doesn't work here
  Copt.opts.Dscale = ones(size(Copt.Dscale));
  while ( (lend-l)*Copt.startdir >=0 && n<Copt.Nmax )      
      %% Polynomial Predictor
      Ulpred = PREDICTORFUN(Ulc(:, 1:n), dUlds(:, 1:n), Ss(1:n), ds/Copt.Dscale(end), Copt.predictororder);

      [Ulc(:, n+1), ~, eflag, nits] = GSOLVE(@(ul) COMBRFUN(func, ul, Ulc(:,n), dUlds(:,n), ds/Copt.Dscale(end), Copt.arclengthparm, Copt.Dscale), Ulpred./Copt.Dscale, Copt.opts);
      ierr = 1;
      ds_ref = ds;      
      while eflag<=0 && ds>Copt.dsmin  % Reduce Step-Size
          ds = ds_ref/2^ierr;
          
          Ulpred = PREDICTORFUN(Ulc(:, 1:n), dUlds(:, 1:n), Ss(1:n), ds/Copt.Dscale(end), Copt.predictororder);
          [Ulc(:, n+1), ~, eflag, nits] = GSOLVE(@(ul) COMBRFUN(func, ul, Ulc(:,n), dUlds(:,n), ds/Copt.Dscale(end), Copt.arclengthparm, Copt.Dscale), Ulpred./Copt.Dscale, Copt.opts);  
          if eflag>0
              break;
          end
          
          ds = ds_ref*2^ierr;
          Ulpred = PREDICTORFUN(Ulc(:, 1:n), dUlds(:, 1:n), Ss(1:n), ds/Copt.Dscale(end), Copt.predictororder);
          [Ulc(:, n+1), ~, eflag, nits] = GSOLVE(@(ul) COMBRFUN(func, ul, Ulc(:,n), dUlds(:,n), ds/Copt.Dscale(end), Copt.arclengthparm, Copt.Dscale), Ulpred./Copt.Dscale, Copt.opts);
          ierr = ierr+1;
      end
      Ulc(:, n+1) = Copt.Dscale.*Ulc(:, n+1);
%       plot(Ulpred(2), Ulpred(1), 'k*');      
%       plot(Ulc(2, 1:n+1), Ulc(1, 1:n+1), '.-'); hold on      
      
      %% Update "Dscale"
      if Copt.DynScale
        Copt.Dscale = max(abs(Ulc(:, n+1)), Copt.Dscale);
        Copt.opts.Dscale = Copt.opts.Dscale;
      end      
      
      %% New Tangent
      [~, dRdU, dRdl] = func(Ulc(:, n+1));
      dRdU = dRdU*diag(Copt.Dscale(1:end-1)); dRdl = dRdl*Copt.Dscale(end);
      z = -dRdU\dRdl;
      al = sign([z; 1]'*dUlds(:, n))/sqrt(1+z'*z);
      dUlds(:, n+1) = [z; 1]*al;

      if eflag<=0
          disp('No convergence');
          break;
      end
      l = Ulc(end, n+1);
      Ss(n+1) = Ss(n) + ds;
      
      %% Reset/Adapt ds
      ds = ds_ref;
      if Copt.stepadapt && prod(dUlds(end, n:n+1))>0
           if nits>Copt.itopt
               ds = max(ds/2, Copt.dsmin);
           else
               ds = min(ds*2, Copt.dsmax);
           end
      end
      
      %% Update n
      n = n+1;
      %% Print
      if Copt.Display
        fprintf('%d %f %f %d: %d\n', n, l, ds, nits, eflag);
      end
      
  end
  
  Ulc = Ulc(:, 1:n);
  dUlds = dUlds(:, 1:n);
  Ss = Ss(1:n);
end

%% Predictor Function
function [Ulpred] = PREDICTORFUN(Uls, dUls, Ss, ds, order)
%%PREDICTORFUN returns the polynomial predictor of appropriate order    
    switch order
        case 0 % Local Predictor
            Ulpred = Uls(:, end);
            Ulpred(end) = Ulpred(end)+ds;
        case 1  % Linear Predictor
            Ulpred = Uls(:, end) + dUls(:, end)*ds;
        otherwise  % Polynomial extrapolator
            order = floor(order/2)*2+1;  % Need to be of odd order
            pts = (order-1)/2+1;
            while pts>length(Ss)
                order = order-2;
                pts = (order-1)/2+1;
            end
            
            Sps = Ss((end-pts+1):end);
            Mx = zeros(order+1, order+1);
            for pi=1:pts
                Mx(pi, :) = Sps(pi).^(0:order);  % Values
                
                Mx(pts+pi, :) = [0 (1:order).*(Sps(pi).^(0:(order-1)))];  % Gradients
            end
            COFS = [(Ss(end)+ds).^(0:order)]*inv(Mx);
            
            Ulpred = [Uls(:, (end-pts+1):end) dUls(:, (end-pts+1):end)]*COFS';
    end
end

%% Combined Residual 
function [R, dRdUl] = COMBRFUN(func, ul, ul0, dulds0, ds, parm, Dscale)
%     [R, dRdU, dRdl] = func(Dscale.*ul);
%     [g, dgdU, dgdl] = ALPARM(ul, ul0./Dscale, dulds0, ds, parm);
%     
%     R = [R; g];
%     dRdUl = [[dRdU dRdl]*diag(Dscale);
%              dgdU dgdl];

    [R, dRdU, dRdl] = func(Dscale.*ul);
    [g, dgdU, dgdl] = ALPARM(Dscale.*ul, ul0, dulds0, ds, parm);
    
    R = [R; g];
    dRdUl = [[dRdU dRdl]*diag(Dscale);
             [dgdU dgdl]*diag(Dscale)];
end

%% Arc Length Parametrization Function
function [g, dgdU, dgdlam] = ALPARM(ul, ul0, dudl0, ds, parm)
    switch parm
        case {'arclength', 'Arclength'}
            g      = (ul-ul0)'*(ul-ul0)-ds^2;
            dgdU   = 2*(ul(1:end-1)-ul0(1:end-1))';
            dgdlam = 2*(ul(end)-ul0(end));
            
%             g      = sign(dudl0(end))*sqrt((ul-ul0)'*(ul-ul0))-ds;
%             dgdU   = 0.5*sign(dudl0(end))/sqrt((ul-ul0)'*(ul-ul0))*dgdU;
%             dgdlam = 0.5*sign(dudl0(end))/sqrt((ul-ul0)'*(ul-ul0))*dgdlam;
        case {'orthogonal', 'Orthogonal'}
            g      = (ul-ul0)'*dudl0-ds;
            dgdU   = dudl0(1:end-1)';
            dgdlam = dudl0(end);
        otherwise
            error('Unknown parametrization');
        end
end