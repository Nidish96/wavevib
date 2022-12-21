function [Ulc, dUlds, Ss, Outs] = PRECOCONT(func, u0, lstart, lend, ds, varargin)
%%PRECONT Uses a bordering update for continuation.
%
%   USAGE:
%       [Uc] = PRECOCONT(func, u0, l0, l1, dl, varargin)
%   INPUTS:
%       func, u0, l0, l1, dl, Copt
%   OUTPUTS:
%       Uc

  %% Default options
  Copt = struct('Nmax', 100, 'dsmax', ds*5, 'dsmin', ds/5, 'stepadapt', 1, ...
                'angopt', pi/6, 'startdir', sign(lend-lstart), 'dlmax', abs(lend-lstart)/10, ...
                'Display', true, 'itopt', 5, ...
                'itDisplay', false, 'lsrch', 0, 'lsit', 10, ...
                'predictororder', 1, ...  % 1, 3, 5, ...
                'arclengthparm', 'orthogonal',...  % 'orthogonal', 'arclength'
                'Dscale', ones(length(u0)+1,1), ...  % Parameter weight
                'CallbackFun', [], ...
                'ITMAX', 10, 'DynScale', 0, 'crit', 14, ...
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
  Copt.opts.crit = Copt.crit;
  
  %% 
  if Copt.predictororder>3
      warning('High order predictors are known to perform badly around turning points')
  end
  
  %% Allocations
  Ss = zeros(Copt.Nmax, 1);
  Ulc = zeros(length(u0)+1, Copt.Nmax);
  dUlds = zeros(length(u0)+1, Copt.Nmax);
  if ~isempty(Copt.CallbackFun)
      Outs = cell(Copt.Nmax,1);
  else
      Outs = [];
  end

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
  
  %% Dynamic Scaling
  if Copt.DynScale
    Copt.Dscale = max(abs(Ulc(:, 1)), Copt.Dscale);
  end  
  Copt.opts.Dscale = Copt.Dscale;  
  
  %% Initial Tangent
  [R, dRdU, dRdl] = func(Ulc(:,1));
  dRdU = dRdU*diag(Copt.Dscale(1:end-1));   dRdl = dRdl*Copt.Dscale(end);
  z  = -dRdU\dRdl;
  al = Copt.startdir/sqrt(1+z'*z);
  dUlds(:, 1) = [z; 1]*al;  % Normalized in scaled-space but expressed in physical space
  dUlds(:, 1) = Copt.Dscale.*dUlds(:,1);
  
  %% BEGIN CONTINUATION
  l = lstart;
  n = 1;
  
  if Copt.Display
      fprintf('N Lam ds Its: Flag\n');
  end
  while ( (lend-l)*Copt.startdir >=0 && n<Copt.Nmax )      
      %% Polynomial Predictor
      Ulpred = PREDICTORFUN(Ulc(:, 1:n), dUlds(:, 1:n), Ss(1:n), ds, Copt.predictororder);
      [Ulc(:, n+1), eflag, dUlds(:, n+1), nits] = CORRECTORFUN(func, Ulc(:,n)./Copt.Dscale, Ulpred./Copt.Dscale, dUlds(:,n)./Copt.Dscale, ds, Copt.arclengthparm, Copt.opts);

    %% Deal with Non-Convergence
    if eflag<=0         
        % 1. Try with Secant Predictor
        Ulpred = PREDICTORFUN(Ulc(:, 1:n), dUlds(:, 1:n), Ss(1:n), ds, -1);  % Secant Predictor
        tmpDer = (Ulpred-Ulc(:,n))/ds;  tmpDer = tmpDer/vecnorm(tmpDer./Copt.Dscale);
        [Ulc(:, n+1), eflag, dUlds(:, n+1), nits] = CORRECTORFUN(func, Ulc(:,n)./Copt.Dscale, Ulpred./Copt.Dscale, tmpDer./Copt.Dscale, ds, Copt.arclengthparm, Copt.opts);

        if eflag<=0
            % 2. Try Natural Parameterization at Predicted Point
            Ulpred = PREDICTORFUN(Ulc(:, 1:n), dUlds(:, 1:n), Ss(1:n), ds, Copt.predictororder);
            Ulc(end,n+1) = Ulpred(end);
            [Ulc(1:end-1,n+1), ~, eflag, ~] = GSOLVE(@(u) func([u; Ulc(end,n+1)]), Ulpred(1:end-1), Copt.opts);  % Refine predictor 
            if eflag>0
                [~, dRdU, dRdl] = func(Ulc(:,n+1));
                dRdU = dRdU*diag(Copt.Dscale(1:end-1));     dRdl = dRdl*Copt.Dscale(end);
                z = -dRdU\dRdl;
                al = sign((Copt.Dscale.*[z; 1])'*dUlds(:, n))/sqrt(1+z'*z);
                dUlds(:, n+1) = [z; 1]*al;
                dUlds(:, n+1) = Copt.Dscale.*dUlds(:, n+1);
            end
            
            if eflag<=0
                % 3. Try refining step size
                ierr = 1;
                ds_ref = ds;
                while eflag<=0 && ds>Copt.dsmin  % Reduce Step-Size
                    ds = ds_ref/2^ierr;

                    Ulpred = PREDICTORFUN(Ulc(:, 1:n), dUlds(:, 1:n), Ss(1:n), ds/Copt.Dscale(end), Copt.predictororder);
                    %           [Ulpred(1:end-1), ~, ~, ~] = GSOLVE(@(u) func([u; Ulpred(end)]), Ulpred(1:end-1), Copt.opts);  % Refine predictor
                    [Ulc(:, n+1), eflag, dUlds(:, n+1), nits] = CORRECTORFUN(func, Ulc(:,n)./Copt.Dscale, Ulpred./Copt.Dscale, dUlds(:,n), ds/Copt.Dscale(end), Copt.arclengthparm, Copt.opts);
                    if eflag>0
                        break;
                    end

                    ds = ds_ref*2^ierr;
                    Ulpred = PREDICTORFUN(Ulc(:, 1:n), dUlds(:, 1:n), Ss(1:n), ds/Copt.Dscale(end), Copt.predictororder);
                    %           [Ulpred(1:end-1), ~, ~, ~] = GSOLVE(@(u) func([u; Ulpred(end)]), Ulpred(1:end-1), Copt.opts);  % Refine predictor
                    [Ulc(:, n+1), eflag, dUlds(:, n+1), nits] = CORRECTORFUN(func, Ulc(:,n), Ulpred./Copt.Dscale, dUlds(:,n), ds/Copt.Dscale(end), Copt.arclengthparm, Copt.opts);          
                    ierr = ierr+1;
                end
            end
        end

        % Exit if nothing has worked so far
        if eflag<=0
            disp('No convergence');
            break;
        end
    end

    if dUlds(end, n+1)*dUlds(end, n)<0 % Gradient Direction has changed - reset gradient to finite difference estimate - this is done to ensure we step through even though there might be a bifurcation
%         keyboard
        dUlds(:, n+1) = (Ulc(:, n+1)-Ulc(:, n));
        dUlds(:, n+1) = dUlds(:, n+1)/vecnorm(dUlds(:, n+1)./Copt.Dscale);
    end      
    
    if ~isempty(Copt.CallbackFun)
        Outs{n+1} = Copt.CallbackFun(Ulc(:, 1:n+1), dUlds(:, 1:n+1), Ulpred, Copt);
        drawnow;
    end
      
      %% Accept Solution Point
      l = Ulc(end, n+1);
      Ss(n+1) = Ss(n) + ds;
      
      %% Update "Dscale"
      if Copt.DynScale
        Copt.Dscale = max(abs(Ulc(:, n+1)), Copt.Dscale);
        Copt.opts.Dscale = Copt.Dscale;
      end
      
      %% Reset/Adapt ds
     if Copt.stepadapt && prod(dUlds(end, n:n+1))>0
          xi = min(max((Copt.itopt/nits)^2, 0.5), 2);
          
          ds = min(max(ds*xi, Copt.dsmin), Copt.dsmax);
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
  if ~isempty(Copt.CallbackFun)
      Outs = Outs{1:n};
  end
end

%% Predictor Function
function [Ulpred] = PREDICTORFUN(Uls, dUls, Ss, ds, order)
%%PREDICTORFUN returns the polynomial predictor of appropriate order    
    switch order
        case -1  % Secant predictor
            if length(Ss)==1
                Ulpred = [Uls(1:end-1, end); Uls(end, end)+ds];
            else
                Ulpred = Uls(:, end) + ds/(Ss(end)-Ss(end-1))*(Uls(:, end)-Uls(:, end-1));
            end
        case 0 % Zero order Predictor
            Ulpred = [Uls(1:end-1, end); Uls(end, end)+ds];
        case 1  % Linear Predictor
            Ulpred = Uls(:, end) + dUls(:, end)*ds;
        otherwise  % Polynomial extrapolator (Hermite): The gradients are not scaled correctly so this has to be fixed before using this
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

%% Corrector Function
function [Ul, eflag, dUl, nits] = CORRECTORFUN(func, Ul0, Ulpred, dUl0, ds, parm, opts)
%%CORRECTORFUN Conducts the corrector iterations

    Dscale = opts.Dscale;
    DU = diag(Dscale(1:end-1));
    Dl = Dscale(end);
    
    Ul = Ulpred;
    %% Flag Function
    flagfun = @(e, e0, r, u) 1*(u<opts.utol) + 2*(e/e0<opts.reletol) + ...
        4*(e<opts.etol) + 8*(r<opts.rtol) + ...
        16*(u<eps | r<eps | abs(e/e0)<eps);

    %% Evaluate Algebraic function & Arc-Length Parametrization
    [R, dRdU, dRdl] = func(Dscale.*Ul);
    [gfun, dgU, dgl] = ALPARM(Ul, Ul0, dUl0, ds, parm);
    % Apply scaling to gradients
    dRdU = dRdU*DU;     dRdl = dRdl*Dl;
    
    %% Bordering Updates Updates
    z12 = dRdU\[-R -dRdl];
    if any(~isfinite(z12(:)))
        z12 = pinv(dRdU)*[-R -dRdl];
    end
    dl = -(gfun+dgU*z12(:,1))/(dgl+dgU*z12(:,2));
    dU = z12(:,1)+z12(:,2)*dl;
    
    %% Error Norms
    r0 = vecnorm([R; gfun]);
    u0 = vecnorm([dU; dl]);
    e0 = abs([R; gfun]'*[dU; dl]);
    eflag = flagfun(e0, e0, r0, u0);
    
	it = 0;
    if opts.Display
        fprintf('ITN, E, E/E0, r, du: eflag\n');
        fprintf('%d, %e, %e, %e, %e: %d\n', 0, e0, 1.0, r0, u0, eflag);
    end
    
    %% Carry Out Corrector
    while (eflag<=opts.crit || it<1) && (it<opts.ITMAX)
        %% Use Line Search if necessary
        if opts.lsrch~=0  % Line-search
            ls_s = 1.0;
            [R, dRdU, dRdl]  = func(Dscale.*(Ul+ls_s*[dU;dl]));
            [gfun, dgU, dgl] = ALPARM(Ul+ls_s*[dU;dl], Ul0, dUl0, ds, parm);
            dRdU = dRdU*DU;     dRdl = dRdl*Dl;

            rp = [R; gfun]'*[R; gfun];
            drdsp = 2*[R; gfun]'*[dRdU*dU+dRdl*dl; dgU*dU+dgl*dl];
            dls = -rp/drdsp;
            if ~isfinite(dls) || ds<=0
                dls = 0;
            end
            for li=1:opts.lsit
                [R, dRdU, dRdl]  = func(Ul+(ls_s+dls)*[DU*dU;Dl*dl]);
                [gfun, dgU, dgl] = ALPARM(Ul+(ls_s+dls)*[dU;dl], Ul0, dUl0, ds, parm);
                dRdU = dRdU*DU;     dRdl = dRdl*Dl;
                
                r = [R; gfun]'*[R; gfun];
                if r>=rp || abs(dls)<1e-10
                    break;
                else
                    ls_s = ls_s+dls;
                end
                drds = 2*[R; gfun]'*[dRdU*dU+dRdl*dl; dgU*dU+dgl*dl];
                dls = -r/drds;
                
                rp = r;
            end
            
            dU = ls_s*dU;
            dl = ls_s*dl;
        end
        [R, dRdU, dRdl]  = func(Dscale.*(Ul+[dU;dl]));
        [gfun, dgU, dgl] = ALPARM(Ul+[dU;dl], Ul0, dUl0, ds, parm);
        dRdU = dRdU*DU;     dRdl = dRdl*Dl;
        
        %% Conduct Update
        Ul = Ul + [dU; dl];
        it = it+1;
        
        %% Bordering Updates
        z12 = dRdU\[-R -dRdl];
        if any(~isfinite(z12(:)))
            z12 = pinv(dRdU)*[-R -dRdl];
        end
        dl = - (gfun+dgU*z12(:,1))/(dgl+dgU*z12(:,2));
        dU = z12(:,1)+z12(:,2)*dl;
        
        %% Termination criteria
        r = vecnorm([R; gfun]);
        u = vecnorm([dU; dl]);
        e = abs([R; gfun]'*[dU; dl]);
        
        eflag = flagfun(e, e0, r, u);
        
        if opts.Display
            fprintf('%d, %e, %e, %e, %e: %d\n', it, e, e/e0, r, u, eflag);
        end
    end
    
    %% Tangent 
    z = -dRdU\dRdl;
    al = sign([z; 1]'*dUl0)/sqrt(1+z'*z);
    dUl = [z; 1]*al;
    nits = it;
    
    %% Scale Back Solution
    Ul = Dscale.*Ul;
    dUl = Dscale.*dUl;
end

%% Combined Residual 
function [R, dRdUl] = COMBRFUN(func, ul, ul0, dudl0, ds, parm, Dscale)    
    [R, dRdU, dRdl] = func(Dscale.*ul);
    [g, dgdU, dgdl] = ALPARM(ul, ul0./Dscale, dudl0./Dscale, ds, parm);
    
    R = [R; g];
    dRdUl = [[dRdU dRdl]*diag(Dscale);
             dgdU dgdl];
end

%% Arc Length Parametrization Function
function [g, dgdU, dgdlam] = ALPARM(ul, ul0, dudl0, ds, parm)
    dudl0 = dudl0/vecnorm(dudl0);

    switch parm
        case {'arclength', 'Arclength'}
            g      = (ul-ul0)'*(ul-ul0)-ds^2;
            dgdU   = 2*(ul(1:end-1)-ul0(1:end-1))';
            dgdlam = 2*(ul(end)-ul0(end));
        case {'orthogonal', 'Orthogonal'}
            g      = dudl0'*(ul-ul0)-ds;
            dgdU   = dudl0(1:end-1)';
            dgdlam = dudl0(end);
        otherwise
            error('Unknown parametrization');
        end
end
