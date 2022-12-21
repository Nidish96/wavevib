function [U, SKIPSTACK, R, J] = CONTINUEf(func, u0, lam0, lam1, ds, varargin)
%CONTINUE Conducts the continuation and solves the system
%
% USAGE:
% INPUTS:
% OUTPUTS:
    % Default options
	Copt = struct('Nmax', 100, 'tole', 1e-6, 'tolr', 1e-6, 'ITMAX', 1000, ...
        'dsmax', ds*5, 'dsmin', ds/5, 'angopt', pi/6, 'startdir', 1,...
        'Display', 1, 'nev', 1, 'adj', 1, 'zsc', 1.0, ...
        'Lvec', ones(length(u0)+1,1), 'Rvec', ones(length(u0)+1,1), ...
        'opt', optimoptions('fsolve', 'SpecifyObjectiveGradient', false, 'Display', 'off',...
        'MaxIterations', 100));
    if nargin==6
        nflds = fieldnames(varargin{1});
        for i=1:length(nflds)
            Copt.(nflds{i}) = varargin{1}.(nflds{i});
        end
    end

    % Allocations
    U = zeros(length(u0)+1, Copt.Nmax);
    SKIPSTACK.U = zeros(length(u0)+1, fix(Copt.Nmax/10));
    if nargout>=3
        R = zeros(length(u0), Copt.Nmax);
        SKIPSTACK.R = zeros(length(u0), fix(Copt.Nmax/10));
    end
    if nargout>=4
        J = zeros(length(u0), length(u0)+1, Copt.Nmax);
        SKIPSTACK.J = zeros(length(u0), length(u0)+1, fix(Copt.Nmax/10));
    end
    if nargout==1
        clear SKIPSTACK
    end
    
    % Correct initial solution
    [u0s, ~, eflag] = fsolve(@(u) FRES(func, [u; lam0], Copt.Rvec(1:end-1), Copt.Lvec(1:end-1)), u0./Copt.Rvec(1:end-1), Copt.opt);
    u0s = Copt.Rvec(1:end-1).*u0s;
    if eflag <= 0
        error('Initial point non-convergent!');
    elseif Copt.Display
        disp('Initial Point Converged')
    end
    U(:, 1) = [u0s; lam0];
    [R0, dRdX0, dRdlam0] = func(U(:, 1));
    if nargout>=3
        R(:, 1) = R0;
    end
    if nargout>=4
        J(:, :, 1) = [dRdX0 dRdlam0];
    end
    
    % BEGIN CONTINUATION
    lam  = lam0;  % current parameter
    lamp = lam0;  % previous parameter
    n    = 1;
    
    % Initial tangent
    z = -dRdX0\dRdlam0;
    z(~isfinite(z)) = 0;
    al = 1.0/sqrt(1+z'*z)*Copt.startdir;
    
    alp = al;
    zp  = z;
    
    duds   = al*[z; 1];
    uguess = U(:, 1) + ds*duds;
    skipnext = 0;
    while ( (lam-lam1)*(lamp-lam1) >= 0 && n<Copt.Nmax )
        [U(:, n+1), Rtmp, eflag, out, Jtmp] = fsolve(@(u) EXRES(func, u, U(:, n), duds, ds, Copt.Rvec, Copt.Lvec), uguess./Copt.Rvec, Copt.opt);
        U(:, n+1) = Copt.Rvec.*U(:, n+1);
        Jtmp = (1./Copt.Lvec).*Jtmp.*(1./Copt.Rvec');
        if eflag<=0
            if ds == Copt.dsmin
                if max(abs(uguess))<eps
                    disp('Diverged!');
                    break;
                else
                    disp('Diverged! Trying with zero initial guess!');
                    uguess = uguess*0;
                    continue;
                end
            else
               disp('Diverged - reducing step-size');
               ds = max(ds/2,Copt.dsmin);
               uguess = U(:, n) + ds*duds;
               continue;
           end
        end
        if nargout>=3
            R(:, n+1)    = (1./Copt.Lvec(1:end-1)).*Rtmp(1:end-1);
        end
        if nargout>=4
            J(:, :, n+1) = Jtmp(1:end-1,:);
        end
        
        % Step
        lamp = lam;
        lam  = U(end, n+1);
        
        % tangent and predictor
        z = -Jtmp(1:end-1, 1:end-1)\Jtmp(1:end-1, end);
        z(~isfinite(z)) = 0;
        al = 1.0/sqrt(1+z'*z) * sign(alp*(1+zp'*z));
        
        % step size adaptation
        theta = acos(alp*al*(1+Copt.zsc^2*zp'*z));
        if Copt.Display
            fprintf('%d %e %e %e %d\n', n+1, U(end,n+1), ds, theta, out.iterations);
        end
        if abs(theta)>abs(Copt.angopt) && abs(ds)>abs(Copt.dsmin) % optimal angle
            % Save current point in skipstack
            if nargout>=2
                SKIPSTACK.Ns(skipnext+1)   = n+1+skipnext;
                SKIPSTACK.U(:, skipnext+1) = U(:, n+1);
                if nargout>=3
                    SKIPSTACK.R(:, skipnext+1) = R(:, n+1);
                end
                if nargout>=4
                    SKIPSTACK.J(:, :, skipnext+1) = J(:, :, n+1);
                end
                skipnext = skipnext+1;
            end
            
            % reset guess and continue
            ds = max(Copt.dsmin, ds/2);
            uguess = U(:, n) + ds*duds;
            continue;
        elseif out.iterations<=50
            ds = min(Copt.dsmax, ds*2);
        end
        
        alp = al;
        zp  = z;
        
        n = n+1;        
        duds   = al*[z; 1];
        uguess = U(:, n) + ds*duds;
    end
    U = U(:, 1:n);
    if nargout>=2
        SKIPSTACK.U = SKIPSTACK.U(:, 1:skipnext);
        if nargout>=3
            R = R(:, 1:n);
            SKIPSTACK.R = SKIPSTACK.R(:, 1:skipnext);
        end
        if nargout>=4
            J = J(:, :, 1:n);
            SKIPSTACK.J = SKIPSTACK.J(:, :, 1:skipnext);
        end
    end
    
    if eflag>0
        disp('Continuation completed successfully');
    end
end

function [Re, Je] = EXRES(func, u, u0, up0, ds, Rvec, Lvec)
%EXRES is the residual function along with the continuation constraint
    [Re, dRdUe, dRdLe] = func(Rvec.*u);    
    Re = Lvec.*[Re; up0'*(Rvec.*u-u0)-ds];
    Je = Lvec.*([dRdUe dRdLe; up0'].*Rvec');
end

function [R, J] = FRES(func, u, Rvec, Lvec)
    [R, J] = func([Rvec; 1].*u);
    R = Lvec.*R;
    J = Lvec.*(J.*Rvec');
end