function [Amat, dAmatdw, dAmatdxi, Fv, dFvdw, dFvdxi, JEV] = WVAMAT(wxi, h, pcs, bcs, joints, Klib, varargin)
%WVAMAT returns the linear "A" matrix and "F" vector for the wave-based
%model. Also returned are the Jacobians w.r.t. w (frequency) and xi
%(parameter); and selector-projector matrices for evaluating the
%nonlinearities.
% Returns in either complex or real-imaginary representation.
%   This routine is for the general quasi-periodic case (with multiple
%   incommensurate frequency components).
%
%   USAGE:
%       [Amat, dAmatdw, dAmatdxi, Fv, dFvdw, dFvdxi, JEV] = WVAMAT(wxi, h, pcs, bcs, joints, Klib);
%           (OR)
%       [Amat, dAmatdw, dAmatdxi, Fv, dFvdw, dFvdxi, JEV] = WVAMAT(wxi, h, pcs, bcs, joints, Klib, 'r');
%   INPUTS:
%       wxi     : (Nc+Npar,1) vector of (frequencies,parameter)
%                   [w1;w2;...;xi1;x2;...]
%       h       : (Nh,Nc) list of harmonics
%       pcs     : (array of structs) wave-based pieces with fields,
%           'coords'
%           'wcomps'
%           'irange'
%           'exci'
%           'excnh'
%           'exccofs'
%       bcs     : (array of structs) boundary conditions with fields,
%           'i', 'pi'
%           'cofs', 'dcofsdw', 'dcofsdxi'
%           'rhs', 'drhsdw', 'drhsdxi'
%           'rih'
%       joints  : (array of structs) joints' informations with fields,
%           'type' (currently only supports 2)
%           'i', 'j'            , 'pi', 'pj'
%           'cofs', 'dcofsdw', 'dcofsdxi'
%           'rhs', 'drhsdw', 'drhsdxi'
%               If nonlinear,
%           'nl',
%           'nldcofs', 'dnldcofsdw', 'dnldcofsdxi'
%           'nlfcofs', 'dnlfcofsdw', 'dnlfcofsdxi'
%           'nld'
%       Klib    : (array of structs) function handles of dispersion
%           relationships in the form of K=f(w,xi) with fields
%           'K', 'dKdw', 'dKdxi'
%           (OPTIONAL)
%       rep     : ('c' or 'r') [Default] 'c' - complex representation
%                                        'r' - real-imaginary
%                                        representation
%   OUTPUTS:
%       Amat    : (Npts*Nwc*Nh, Npts*Nwc*Nh) OR (Npts*Nwc*Nhc, Npts*Nwc*Nhc)
%       dAmatdw : (Npts*Nwc*Nh, Npts*Nwc*Nh) OR (Npts*Nwc*Nhc, Npts*Nwc*Nhc)
%       dAmatdxi: (Npts*Nwc*Nh, Npts*Nwc*Nh) OR (Npts*Nwc*Nhc, Npts*Nwc*Nhc)
%       Fv      : (Npts*Nwc*Nh, 1) OR (Npts*Nwc*Nhc, 1)
%       dFvdw   : (Npts*Nwc*Nh, 1) OR (Npts*Nwc*Nhc, 1)
%       dFvdxi  : (Npts*Nwc*Nh, 1) OR (Npts*Nwc*Nhc, 1)
%       JEV     : (array of structs) with fields
%           'Lj', 'dLjdw', 'dLjdxi'
%           'Gj', 'dGjdw', 'dGjdxi'


    Npts = pcs(end).irange(end);
    Nwc  = size(pcs(1).wcomps,1);  % Number of wave components
    Nh   = size(h,1);  % Number of harmonic terms
    Nc   = size(h,2);  % Number of incommensurate frequency components
    Npar = length(wxi)-Nc;  % Number of free parameters (xis)

    Amat = zeros(Npts*Nwc*Nh);
    dAmatdw = zeros(Npts*Nwc*Nh,Npts*Nwc*Nh,Nc);
    dAmatdxi = zeros(Npts*Nwc*Nh,Npts*Nwc*Nh,Npar);
    Fv = zeros(Npts*Nwc*Nh,1);
    dFvdw = zeros(Npts*Nwc*Nh,Nc);
    dFvdxi = zeros(Npts*Nwc*Nh,Npar);

    % Setup Selector-Projectors for Nonlinearities
    nlis = find(arrayfun(@(j) ~isempty(j.nl), joints));
    nnl = length(nlis);
    Lj = cell(nnl,1);
    dLjdw = cell(nnl,1);
    dLjdxi = cell(nnl,1);
    Gj = cell(nnl,1);
    dGjdw = cell(nnl,1);
    dGjdxi = cell(nnl,1);
    for n=1:length(nlis)
        k = nlis(n);
        
        Lj{n} = zeros(joints(k).nld*Nh, Npts*Nwc*Nh);
        dLjdw{n} = zeros(joints(k).nld*Nh, Npts*Nwc*Nh,Nc);
        dLjdxi{n} = zeros(joints(k).nld*Nh, Npts*Nwc*Nh,Npar);

        Gj{n} = zeros(Npts*Nwc*Nh, joints(k).nld*Nh);
        dGjdw{n} = zeros(Npts*Nwc*Nh, joints(k).nld*Nh,Nc);
        dGjdxi{n} = zeros(Npts*Nwc*Nh, joints(k).nld*Nh,Npar);
    end

    for ih=1:Nh
        hstart = (ih-1)*(Npts*Nwc);
        
        k = 0;
        % 1. Propagation in Pieces
        for i=1:length(pcs)
            % Wavenumbers of the different components
            Ks = arrayfun(@(k) k.K(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))), Klib);  %(nk,1)
            dKdws = cell2mat(arrayfun(@(k) k.dKdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*h(ih, :), Klib, 'UniformOutput', false));  %(nk,Nc)
            dKdxis = arrayfun(@(k) k.dKdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))), Klib);  %(nk,Npar)

            K = pcs(i).wcomps(:,1).*Ks(pcs(i).wcomps(:,2));  %(Nwc,1)
            dKdw = pcs(i).wcomps(:,1).*dKdws(pcs(i).wcomps(:,2),:);  %(Nwc,Nc)
            dKdxi = pcs(i).wcomps(:,1).*dKdxis(pcs(i).wcomps(:,2),:);  %(Nwc,Npar)

            for n=1:pcs(i).N-1
                k = k+1;
                si = hstart + (k-1)*Nwc+1;
                se = hstart + k*Nwc;
    
                qi = hstart + ((pcs(i).irange(1)-1)+n-1)*Nwc+1;
                qe = hstart + ((pcs(i).irange(1)-1)+n+1)*Nwc;
                
                dx = abs(diff(pcs(i).U(n:n+1))*pcs(i).S);
    
                Amat(si:se,qi:qe) = [diag(exp(K*dx)) -eye(Nwc)];
                dAmatdw(si:se,qi:(qe-Nwc),:) = reshape(cell2mat(arrayfun(@(a) diag(exp(K*dx)*dx.*dKdw(:,a)), 1:Nc, 'UniformOutput', false)), Nwc,Nwc,Nc);
                dAmatdxi(si:se,qi:(qe-Nwc),:) = reshape(cell2mat(arrayfun(@(a) diag(exp(K*dx)*dx.*dKdxi(:,a)), 1:Npar, 'UniformOutput', false)), Nwc,Nwc,Npar);
                
                % Check if this is an excitation Point
                ei = find(pcs(i).exci==n);
                if dx==0 && ~isempty(ei)
                    if length(ei)>1
                        error('More than 1 excitation location detected. Check pcs.');
                    end
    
                    if all(pcs(i).excnh(ei,:)==h(ih,:))
                        Fv(si:se) = pcs(i).exccofs{ei}{1}(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  %(Nwc,1)
                        dFvdw(si:se,:) = pcs(i).exccofs{ei}{2}(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))).*h(ih,:);  %(Nwc,Nc)
                        dFvdxi(si:se,:) = pcs(i).exccofs{ei}{3}(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  %(Nwc,Npar)
                    end
                end
            end
        end

        % 2. Boundary Conditions
        ktn = hstart + (Npts-length(pcs))*Nwc;
        for n=1:length(bcs)
            inds = hstart + (bcs(n).i-1)*Nwc+(1:Nwc);
    
            Amat(ktn+n, inds) = bcs(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  % (1,1)
            dAmatdw(ktn+n, inds, :) = bcs(n).dcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))).*permute(h(ih,:),[1 3 2]);  %(1,Nc)
            dAmatdxi(ktn+n, inds, :) = bcs(n).dcofsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  %(1,Npar)

            % RHS
            if ~isempty(bcs(n).rih) && all(bcs(n).rih==h(ih,:))
                Fv(ktn+n) = Fv(ktn+n)+bcs(n).rhs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));
                dFvdw(ktn+n, :) = dFvdw(ktn+n)+bcs(n).rhs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))).*h(ih,:);
                dFvdxi(ktn+n, :) = dFvdxi(ktn+n)+bcs(n).rhs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));
            end
        end
    
        % 3. Joints
        ktn = hstart + (Npts-length(pcs))*Nwc+length(bcs);
        for n=1:length(joints)
            switch joints(n).type
                case {1, 2}  % Binary Connection
                    inds = hstart + [(joints(n).i-1)*Nwc+(1:Nwc) (joints(n).j-1)*Nwc+(1:Nwc)];
    
                    Amat(ktn+(n-1)*Nwc+(1:Nwc), inds) = joints(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  % (Nwc,2*Nwc)
                    dAmatdw(ktn+(n-1)*Nwc+(1:Nwc), inds, :) = reshape(cell2mat(arrayfun(@(a) joints(n).dcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*h(ih,a), ...
                        1:Nc, 'UniformOutput', false)), Nwc, 2*Nwc, Nc);  % (Nwc,2*Nwc,Nc)
                    dAmatdxi(ktn+(n-1)*Nwc+(1:Nwc), inds, :) = joints(n).dcofsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  %(Nwc,2*Nwc,Npar)
                otherwise
                    error('Needs to be implemented still.');
            end
        end
        % 3a. Nonlinear Joints
        for n=1:nnl
            k = nlis(n);
            nld = joints(k).nld;  % Nonlinear Displacements (& forces)
    
            switch joints(k).type
                case 2
                    % Choosing NL displacement
                    inds = hstart + [(joints(k).i-1)*Nwc+(1:Nwc) (joints(k).j-1)*Nwc+(1:Nwc)];
    
                    Lj{n}((ih-1)*nld+(1:nld), inds) = joints(k).nldcofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  % (nld,Nwc)
                    dLjdw{n}((ih-1)*nld+(1:nld), inds,:) = reshape(cell2mat(arrayfun(@(a) joints(k).dnldcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*h(ih,a), ...
                        1:Nc, 'UniformOutput', false)), nld,2*Nwc,Nc);  %(nld,Nwc,Nc)
                    dLjdxi{n}((ih-1)*nld+(1:nld), inds,:) = joints(k).dnldcofsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  %(nld,Nwc,Npar)
    
                    % Putting NL force
                    kinds = ktn+(k-1)*Nwc+(1:Nwc);
                    
                    Gj{n}(kinds, (ih-1)*nld+(1:nld)) = joints(k).nlfcofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  %(Nwc,nld)
                    dGjdw{n}(kinds, (ih-1)*nld+(1:nld),:) = reshape(cell2mat(arrayfun(@(a) joints(k).dnlfcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*h(ih,a), ...
                        1:Nc, 'UniformOutput', false)), Nwc,nld,Nc);  %(Nwc,nld,Nc)
                    dGjdxi{n}(kinds, (ih-1)*nld+(1:nld),:) = joints(k).dnlfcofsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  %(Nwc,nld,Nc)
                otherwise
                    error('Needs to be implemented still.');
            end
            if ~isempty(joints(k).rih) && all(joints(k).rih==h(ih,:))
                Fv(kinds) = Fv(kinds) + joints(k).rhs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  %(Nwc,1)
                dFvdw(kinds,:) = dFvdw(kinds,:) + joints(k).drhsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))).*h(ih,:);  %(Nwc,Nc)
                dFvdxi(kinds,:) = dFvdxi(kinds,:) + joints(k).drhsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));  % (Nwc,Npar)
            end
        end
    end
    JEV = struct('Lj', Lj, 'dLjdw', dLjdw, 'dLjdxi', dLjdxi, ...
        'Gj', Gj, 'dGjdw', dGjdw, 'dGjdxi', dGjdxi);

    %% Convert to Fully Real Representation
    if length(varargin)>=1 && varargin{1}=='r'
        Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
        [zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);

        Amatc = Amat;
        dAmatdwc = dAmatdw;
        dAmatdxic = dAmatdxi;
    
        Amat = zeros(Npts*Nwc*Nhc);
        dAmatdw = zeros(Npts*Nwc*Nhc, Npts*Nwc*Nhc, Nc);
        dAmatdxi = zeros(Npts*Nwc*Nhc, Npts*Nwc*Nhc, Npar);
    
        Amat([rinds iinds], [rinds iinds]) = [real([Amatc(hinds,hinds) 1j*Amatc(hinds,hinds)]);
            imag([Amatc(hinds,hinds) 1j*Amatc(hinds,hinds)])];
        Amat(rinds0, rinds0) = real(Amatc(zinds,zinds));
        dAmatdw([rinds iinds], [rinds iinds], :) = [real([dAmatdwc(hinds,hinds, :) 1j*dAmatdwc(hinds,hinds, :)]);
            imag([dAmatdwc(hinds,hinds, :) 1j*dAmatdwc(hinds,hinds, :)])];
        dAmatdw(rinds0, rinds0, :) = real(dAmatdwc(zinds,zinds, :));
        dAmatdxi([rinds iinds], [rinds iinds], :) = [real([dAmatdxic(hinds,hinds, :) 1j*dAmatdxic(hinds,hinds, :)]);
            imag([dAmatdxic(hinds,hinds, :) 1j*dAmatdxic(hinds,hinds, :)])];
        dAmatdxi(rinds0, rinds0, :) = real(dAmatdxic(zinds,zinds, :));
    
        Fvc = Fv;
        dFvdwc = dFvdw;
        dFvdxic = dFvdxi;

        Fv = zeros(Npts*Nwc*Nhc,1);
        dFvdw = zeros(Npts*Nwc*Nhc,Nc);
        dFvdxi = zeros(Npts*Nwc*Nhc,Npar);
        Fv([rinds0 rinds iinds]) = [real(Fvc(zinds)); real(Fvc(hinds)); imag(Fvc(hinds))];
        dFvdw([rinds0 rinds iinds], :) = [real(dFvdwc(zinds, :)); real(dFvdwc(hinds, :)); imag(dFvdwc(hinds, :))];
        dFvdxi([rinds0 rinds iinds], :) = [real(dFvdxic(zinds, :)); real(dFvdxic(hinds, :)); imag(dFvdxic(hinds, :))];

        % Nonlinear Selector-Projectors
        JEVc = JEV;
        JEV = struct('Lj', cell(nnl,1), ...
            'dLjdw', cell(nnl,1), ...
            'dLjdxi', cell(nnl,1), ...
            'Gj', cell(nnl,1), ...
            'dGjdw', cell(nnl,1), ...
            'dGjdxi', cell(nnl,1));

        for n=1:nnl
            k = nlis(n);
            nld = joints(k).nld;

            [nlzinds,nlhinds,nlrinds0,nlrinds,nliinds] = HINDS(nld, h);
            JEV(n).Lj = zeros(nld*Nhc, Npts*Nwc*Nhc);
            JEV(n).dLjdw = zeros(nld*Nhc, Npts*Nwc*Nhc, Nc);
            JEV(n).dLjdxi = zeros(nld*Nhc, Npts*Nwc*Nhc, Npar);
            
            JEV(n).Gj = zeros(Npts*Nwc*Nhc, nld*Nhc);
            JEV(n).dGjdw = zeros(Npts*Nwc*Nhc, nld*Nhc, Nc);
            JEV(n).dGjdxi = zeros(Npts*Nwc*Nhc, nld*Nhc, Npar);
           
            tmp = [JEVc(n).Lj(nlhinds, hinds) 1j*JEVc(n).Lj(nlhinds, hinds)];
            JEV(n).Lj(nlrinds0, rinds0) = JEVc(n).Lj(nlzinds, zinds);
            JEV(n).Lj([nlrinds nliinds], [rinds iinds]) = [real(tmp);imag(tmp)];
            tmp = [JEVc(n).dLjdw(nlhinds, hinds, :) 1j*JEVc(n).dLjdw(nlhinds, hinds, :)];
            JEV(n).dLjdw(nlrinds0, rinds0, :) = JEVc(n).dLjdw(nlzinds, zinds, :);
            JEV(n).dLjdw([nlrinds nliinds], [rinds iinds], :) = [real(tmp);imag(tmp)];
            tmp = [JEVc(n).dLjdxi(nlhinds, hinds, :) 1j*JEVc(n).dLjdxi(nlhinds, hinds, :)];
            JEV(n).dLjdxi(nlrinds0, rinds0, :) = JEVc(n).dLjdxi(nlzinds, zinds, :);
            JEV(n).dLjdxi([nlrinds nliinds], [rinds iinds], :) = [real(tmp);imag(tmp)];

            tmp = [JEVc(n).Gj(hinds,nlhinds) 1j*JEVc(n).Gj(hinds,nlhinds)];
            JEV(n).Gj(rinds0,nlrinds0) = JEVc(n).Gj(zinds,nlzinds);
            JEV(n).Gj([rinds iinds],[nlrinds nliinds]) = [real(tmp);imag(tmp)];
            tmp = [JEVc(n).dGjdw(hinds,nlhinds, :) 1j*JEVc(n).dGjdw(hinds,nlhinds, :)];
            JEV(n).dGjdw(rinds0,nlrinds0, :) = JEVc(n).dGjdw(zinds,nlzinds, :);
            JEV(n).dGjdw([rinds iinds],[nlrinds nliinds], :) = [real(tmp);imag(tmp)];
            tmp = [JEVc(n).dGjdxi(hinds,nlhinds, :) 1j*JEVc(n).dGjdxi(hinds,nlhinds, :)];
            JEV(n).dGjdxi(rinds0,nlrinds0, :) = JEVc(n).dGjdxi(zinds,nlzinds, :);
            JEV(n).dGjdxi([rinds iinds],[nlrinds nliinds], :) = [real(tmp);imag(tmp)];
        end
    end
end
