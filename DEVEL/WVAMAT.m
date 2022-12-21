function [Amat, dAmatdw, dAmatdxi, Fv, dFvdw, dFvdxi, JEV] = WVAMAT(wxi, h, pcs, bcs, joints, Klib, wcomps, varargin)
%WVAMAT
    allcoords = cell2mat({pcs.coords}');
    Npts = size(allcoords,1);
    Nwc  = size(wcomps,1);  % Number of wave components
    Nh   = length(h);

    Amat = zeros(Npts*Nwc*Nh);
    dAmatdw = zeros(Npts*Nwc*Nh);
    dAmatdxi = zeros(Npts*Nwc*Nh);
    Fv = zeros(Npts*Nwc*Nh,1);
    dFvdw = zeros(Npts*Nwc*Nh,1);
    dFvdxi = zeros(Npts*Nwc*Nh,1);

    ih = 1;
    % Setup Selector-Projectors for Nonlinearities
    nlis = find(arrayfun(@(j) ~isempty(j.nl), joints));
    nnl = length(nlis);
    Lj = zeros(nnl*Nh, Npts*Nwc*Nh);
    dLjdw = zeros(nnl*Nh, Npts*Nwc*Nh);
    dLjdxi = zeros(nnl*Nh, Npts*Nwc*Nh);
    Gj = zeros(Npts*Nwc*Nh, nnl*Nh);
    dGjdw = zeros(Npts*Nwc*Nh, nnl*Nh);
    dGjdxi = zeros(Npts*Nwc*Nh, nnl*Nh);

    for ih=1:Nh
        hstart = (ih-1)*(Npts*Nwc);
        
        k = 0;
        % Wavenumbers of the different components
        K    = wcomps(:,1).*cellfun(@(c) c(h(ih)*wxi(1),wxi(2)), {Klib(wcomps(:,2)).K}.');
        dKdw = wcomps(:,1).*cellfun(@(c) c(h(ih)*wxi(1),wxi(2)), {Klib(wcomps(:,2)).dKdw}.')*h(ih);
        dKdxi = wcomps(:,1).*cellfun(@(c) c(h(ih)*wxi(1),wxi(2)), {Klib(wcomps(:,2)).dKdxi}.');

        % 1. Propagation in Pieces
        for i=1:length(pcs)
            for n=1:pcs(i).N-1
                k = k+1;
                si = hstart + (k-1)*Nwc+1;
                se = hstart + k*Nwc;
    
                qi = hstart + ((pcs(i).irange(1)-1)+n-1)*Nwc+1;
                qe = hstart + ((pcs(i).irange(1)-1)+n+1)*Nwc;
                
                dx = abs(diff(pcs(i).U(n:n+1))*pcs(i).S);
    
                Amat(si:se,qi:qe) = [diag(exp(K*dx)) -eye(Nwc)];
                dAmatdw(si:se,qi:qe) = [diag(exp(K*dx).*dKdw*dx) zeros(Nwc)];
                dAmatdxi(si:se,qi:qe) = [diag(exp(K*dx).*dKdxi*dx) zeros(Nwc)];
                
                % Check if this is an excitation Point
                if dx==0 && ismember(n,pcs(i).exci)
                    ei = find(pcs(i).exci==n);
                    if length(ei)>1
                        error('More than 1 excitation location detected. Check pcs.');
                    end
    
                    if pcs(i).excnh(ei)==h(ih)
                        Fv(si:se) = pcs(i).exccofs{ei}{1}(h(ih)*wxi(1), wxi(2));
                        dFvdw(si:se) = pcs(i).exccofs{ei}{2}(h(ih)*wxi(1), wxi(2))*h(ih);
                        dFvdxi(si:se) = pcs(i).exccofs{ei}{3}(h(ih)*wxi(1), wxi(2));
                    end
                end
            end
        end

        % 2. Boundary Conditions
        ktn = hstart + (Npts-length(pcs))*Nwc;
        for n=1:length(bcs)
            inds = hstart + (bcs(n).i-1)*Nwc+(1:Nwc);
    
            Amat(ktn+n, inds) = bcs(n).cofs(h(ih)*wxi(1), wxi(2));
            dAmatdw(ktn+n, inds) = bcs(n).dcofsdw(h(ih)*wxi(1), wxi(2))*h(ih);
            dAmatdxi(ktn+n, inds) = bcs(n).dcofsdxi(h(ih)*wxi(1), wxi(2));
        end
    
        % 3. Joints
        ktn = hstart + (Npts-length(pcs))*Nwc+length(bcs);
        for n=1:length(joints)
            switch joints(n).type
                case {1, 2}  % Binary Connection
                    inds = hstart + [(joints(n).i-1)*Nwc+(1:Nwc) (joints(n).j-1)*Nwc+(1:Nwc)];
    
                    Amat(ktn+(n-1)*Nwc+(1:Nwc), inds) = joints(n).cofs(h(ih)*wxi(1), wxi(2));
                    dAmatdw(ktn+(n-1)*Nwc+(1:Nwc), inds) = joints(n).dcofsdw(h(ih)*wxi(1), wxi(2))*h(ih);
                    dAmatdw(ktn+(n-1)*Nwc+(1:Nwc), inds) = joints(n).dcofsdxi(h(ih)*wxi(1), wxi(2));
                otherwise
                    error('Needs to be implemented still.');
            end
        end
        % 3a. Nonlinear Joints
        hstnl = (ih-1)*nnl;
        for n=1:nnl
            k = nlis(n);
    
            switch joints(k).type
                case 2  % Relative Displacement is joint
                    % Choosing NL displacement
                    inds = hstart + [(joints(k).i-1)*Nwc+(1:Nwc) (joints(k).j-1)*Nwc+(1:Nwc)];
    
                    Lj(hstnl+n, inds) = joints(k).nldcofs(h(ih)*wxi(1), wxi(2));
                    dLjdw(hstnl+n, inds) = joints(k).dnldcofsdw(h(ih)*wxi(1), wxi(2))*h(ih);
                    dLjdxi(hstnl+n, inds) = joints(k).dnldcofsdxi(h(ih)*wxi(1), wxi(2));
    
                    % Putting NL force
                    kinds = ktn+(k-1)*Nwc+(1:Nwc);
                    
                    Gj(kinds, hstnl+n) = joints(k).nlfcofs(h(ih)*wxi(1), wxi(2));
                    dGjdw(kinds, hstnl+n) = joints(k).dnlfcofsdw(h(ih)*wxi(1), wxi(2))*h(ih);
                    dGjdxi(kinds, hstnl+n) = joints(k).dnlfcofsdxi(h(ih)*wxi(1), wxi(2));
                otherwise
                    error('Needs to be implemented still.');
            end
        end
    end
    JEV = struct('Lj', Lj, 'dLjdw', dLjdw, 'dLjdxi', dLjdxi, ...
        'Gj', Gj, 'dGjdw', dGjdw, 'dGjdxi', dGjdxi);

    %% Convert to Fully Real Representation
    if length(varargin)>=1 && varargin{1}=='r'
        Nhc = sum((h==0)+2*(h~=0));
        [zinds,hinds,rinds0,rinds,iinds] = HINDS(Npts*Nwc, h);
        [nlzinds,nlhinds,nlrinds0,nlrinds,nliinds] = HINDS(nnl, h);

        Amatc = Amat;
        dAmatdwc = dAmatdw;
        dAmatdxic = dAmatdxi;
    
        Amat = zeros(Npts*Nwc*Nhc);
        dAmatdw = zeros(Npts*Nwc*Nhc);
        dAmatdxi = zeros(Npts*Nwc*Nhc);
    
        Amat([rinds iinds], [rinds iinds]) = [real([Amatc(hinds,hinds) 1j*Amatc(hinds,hinds)]);
            imag([Amatc(hinds,hinds) 1j*Amatc(hinds,hinds)])];
        Amat(rinds0, rinds0) = real(Amatc(zinds,zinds));
        dAmatdw([rinds iinds], [rinds iinds]) = [real([dAmatdwc(hinds,hinds) 1j*dAmatdwc(hinds,hinds)]);
            imag([dAmatdwc(hinds,hinds) 1j*dAmatdwc(hinds,hinds)])];
        dAmatdw(rinds0, rinds0) = real(dAmatdwc(zinds,zinds));
        dAmatdxi([rinds iinds], [rinds iinds]) = [real([dAmatdxic(hinds,hinds) 1j*dAmatdxic(hinds,hinds)]);
            imag([dAmatdxic(hinds,hinds) 1j*dAmatdxic(hinds,hinds)])];
        dAmatdxi(rinds0, rinds0) = real(dAmatdxic(zinds,zinds));
    
        Fvc = Fv;
        dFvdwc = dFvdw;
        dFvdxic = dFvdxi;

        Fv = zeros(Npts*Nwc*Nhc,1);
        dFvdw = zeros(Npts*Nwc*Nhc,1);
        dFvdxi = zeros(Npts*Nwc*Nhc,1);
        Fv([rinds0 rinds iinds]) = [real(Fvc(zinds)); real(Fvc(hinds)); imag(Fvc(hinds))];
        dFvdw([rinds0 rinds iinds]) = [real(dFvdwc(zinds)); real(dFvdwc(hinds)); imag(dFvdwc(hinds))];
        dFvdxi([rinds0 rinds iinds]) = [real(dFvdxic(zinds)); real(dFvdxic(hinds)); imag(dFvdxic(hinds))];

        % Nonlinear Selector-Projectors
        JEVc = JEV;
        JEV = struct('Lj', zeros(nnl*Nhc, Npts*Nwc*Nhc), ...
            'dLjdw', zeros(nnl*Nhc, Npts*Nwc*Nhc), ...
            'dLjdxi', zeros(nnl*Nhc, Npts*Nwc*Nhc), ...
            'Gj', zeros(Npts*Nwc*Nhc, nnl*Nhc), ...
            'dGjdw', zeros(Npts*Nwc*Nhc, nnl*Nhc), ...
            'dGjdxi', zeros(Npts*Nwc*Nhc, nnl*Nhc));
        JEV.Lj([nlrinds0 nlrinds nliinds], [rinds0 rinds iinds]) = blkdiag(JEVc.Lj(nlzinds, zinds),[real([JEVc.Lj(nlhinds, hinds) 1j*JEVc.Lj(nlhinds, hinds)]);imag([JEVc.Lj(nlhinds, hinds) 1j*JEVc.Lj(nlhinds, hinds)])]);
        JEV.dLjdw([nlrinds0 nlrinds nliinds], [rinds0 rinds iinds]) = blkdiag(JEVc.dLjdw(nlzinds, zinds),[real([JEVc.dLjdw(nlhinds, hinds) 1j*JEVc.dLjdw(nlhinds, hinds)]);imag([JEVc.dLjdw(nlhinds, hinds) 1j*JEVc.dLjdw(nlhinds, hinds)])]);
        JEV.dLjdxi([nlrinds0 nlrinds nliinds], [rinds0 rinds iinds]) = blkdiag(JEVc.dLjdxi(nlzinds, zinds),[real([JEVc.dLjdxi(nlhinds, hinds) 1j*JEVc.dLjdxi(nlhinds, hinds)]);imag([JEVc.dLjdxi(nlhinds, hinds) 1j*JEVc.dLjdxi(nlhinds, hinds)])]);

        JEV.Gj([rinds0 rinds iinds],[nlrinds0 nlrinds nliinds]) = blkdiag(JEVc.Gj(zinds,nlzinds),[real([JEVc.Gj(hinds,nlhinds) 1j*JEVc.Gj(hinds,nlhinds)]);imag([JEVc.Gj(hinds,nlhinds) 1j*JEVc.Gj(hinds,nlhinds)])]);
        JEV.dGjdw([rinds0 rinds iinds],[nlrinds0 nlrinds nliinds]) = blkdiag(JEVc.dGjdw(zinds,nlzinds),[real([JEVc.dGjdw(hinds,nlhinds) 1j*JEVc.dGjdw(hinds,nlhinds)]);imag([JEVc.dGjdw(hinds,nlhinds) 1j*JEVc.dGjdw(hinds,nlhinds)])]);
        JEV.dGjdxi([rinds0 rinds iinds],[nlrinds0 nlrinds nliinds]) = blkdiag(JEVc.dGjdxi(zinds,nlzinds),[real([JEVc.dGjdxi(hinds,nlhinds) 1j*JEVc.dGjdxi(hinds,nlhinds)]);imag([JEVc.dGjdxi(hinds,nlhinds) 1j*JEVc.dGjdxi(hinds,nlhinds)])]);
    end
end