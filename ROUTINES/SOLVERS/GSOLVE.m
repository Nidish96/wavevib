function [U, R, eflag, it, dRc, reu] = GSOLVE(func, U0, varargin)
%GSOLVE Uses Gradient-based iterations to solve
%
% USAGE:
%  [U, R, eflag, it, jc] = GSOLVE(func, U0, opts);
% INPUTS:
%  func		: Function handle [F, dFdU] = func(u);
%  U0		: Initial Guess
%  opts		: Options structure with,
% 	reletol (float)
% 	etol    (float)
% 	utol    (float)
% 	rtol    (float)
% 	Display (boolean)
% OUTPUTS:
%  U		:
%  R		:
%  eflag	:
%  it		:
%  Jc		:

  %% Default Parameters
  opts = struct('reletol', 1e-12, 'rtol', 1e-12, 'etol', 1e-12, 'utol', ...
                1e-12, 'Display', false, 'Dscale', ones(size(U0)), ...
	       'ITMAX', 10, ...
           'lsrch', 0, 'lsmin', -10, 'lsmax', 10, 'lsit', 10, ...  % Line Search Parameters
           'crit', 14);  % Termination Criterion - See flagfun below
  if nargin==3
      nflds = fieldnames(varargin{1});
      for i = 1:length(nflds)
	opts.(nflds{i}) = varargin{1}.(nflds{i});
      end
  end
  opts.Dscale = opts.Dscale(:);

  %% Starting Point Residual, Jacobian, and error metrics
  Nu = length(U0);
  U0 = U0./opts.Dscale(1:Nu);
  
  [R0, dR0] = func(opts.Dscale(1:Nu).*U0);
  dU0 = -(dR0*diag(opts.Dscale(1:Nu)))\R0;
  if any(~isfinite(dU0))
      dU0 = -pinv(dR0*diag(opts.Dscale(1:Nu)))*R0;
  end
  e0  = abs(R0'*dU0);
  
  if (e0 < eps)
    e0 = 1.0;
    R0 = ones(Nu, 1);
  end
  r0 = sqrt(R0'*R0);
  u0 = sqrt(dU0'*dU0);
  
  %% Termination Flag
    flagfun = @(e, e0, r, u) 1*(u<opts.utol) + 2*(e/e0<opts.reletol) + ...
        4*(e<opts.etol) + 8*(r<opts.rtol) +...
        16*(u<eps | r<eps | abs(e/e0)<eps);
    % Last condition is just to ensure we're not at machine precision in
    % the updates
    
  %% Initializing for iterations  
  r   = r0;
  e   = e0;
  u   = u0;

  reu = [r e u];
  
  dRc = dR0;
  R   = R0;
  U   = U0;
  dU  = dU0;
  it  = 0;

  eflag = flagfun(e, e0, r, u); 
  if opts.Display
    fprintf('ITN, E, E/E0, r, du\n%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, eflag);
  end
  eflagp = 0;

  %% Start Iterations  
  while (eflag<=opts.crit || it<1) && it<=opts.ITMAX
    %% Use Linear Search if necessary
    if opts.lsrch~=0        
        ls_s = 1;
        [Rp, dRp] = func(opts.Dscale(1:Nu).*(U+dU));
        rp = Rp'*Rp;
        drdsp = 2*Rp'*(dRp*diag(opts.Dscale(1:Nu)))*dU;
        ds = -drdsp\rp;
        if ~isfinite(ds)
            ds=0;
        end
        for li=1:opts.lsit
            [R, dRc] = func(opts.Dscale(1:Nu).*(U+(ls_s+ds)*dU));
            r = R'*R;
            
            if r>=rp || ds<1e-10
                break;
            else
                ls_s = ls_s+ds;
            end
            drds = 2*R'*(dRc*diag(opts.Dscale(1:Nu)))*dU;
            ds = -drds\r;
            if ~isfinite(ds)
                ds = 0;
            end

            rp = r;
            Rp = R;
            dRp = dRc;
        end
        dU = ls_s*dU;
        R = Rp;
        dRc = dRp;
    else
        [R, dRc] = func(opts.Dscale(1:Nu).*(U+dU));
    end
      
	%% Gradient Update    
	U  = U+dU;
    it = it+1;
    
    dU = -(dRc*diag(opts.Dscale(1:Nu)))\R;
    if any(~isfinite(dU))  % Ill-conditioned update: use pinv instead
      dU = -pinv(dRc*diag(opts.Dscale(1:Nu)))*R0;
    end
    
    r = sqrt(R'*R);
    e = abs(R'*dU);
    u = sqrt(dU'*dU);
    
    reu = [reu; [r e u]];

    %% Check Termination
    eflag = flagfun(e, e0, r, u); 
    if opts.Display
      fprintf('%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, eflag);
    end
    
    if it>opts.ITMAX
      eflag = 0;
      break;
    end
  end
  
  % Rescale Solution
  U = opts.Dscale(1:Nu).*U;

  if eflag==0
    disp('No Convergence : Returning')
    return;
  end

  if opts.Display
    disp('Successful Campaign!')
  end
end
