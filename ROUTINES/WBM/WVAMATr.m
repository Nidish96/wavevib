function [Amat, dAmatdw, dAmatdxi, Fv, dFvdw, dFvdxi, JEV] = WVAMATr(wxi, h, pcs, bcs, joints, Klib, varargin)
%WVAMATr returns the linear "A" matrix and "F" vector for the wave-based
%model using reduced parameterization. Also returned are the Jacobians w.r.t.
%w (frequency) and xi (parameter) and selector-projector matrices for
%evaluating the nonlinearities.
% Returns in either complex or real-imaginary representation.
%
%   This parameterization is used for the "smallest-sized" WBMs. Each
%   wave-based piece is represented using exactly one point. This is
%   achieved by solving the linear transmission conditions within each
%   piece. Computationally, this scales better for most problems and has to
%   be preferred over the other version. 
%       The "reduced parameterization" routines are suffixed with "r".
%
%   USAGE:
%       [Amat, dAmatdw, dAmatdxi, Fv, dFvdw, dFvdxi, JEV] = WVAMATr(wxi, h, pcs, bcs, joints, Klib);
%           (OR)
%       [Amat, dAmatdw, dAmatdxi, Fv, dFvdw, dFvdxi, JEV] = WVAMATr(wxi, h, pcs, bcs, joints, Klib, 'r');
%   INPUTS:
%       wxi     : (2,1) vector of (frequency,parameter)
%       h       : (1,Nh) list of harmonics
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
%       Amat    : (Npcs*Nwc*Nh, Npcs*Nwc*Nh) OR (Npcs*Nwc*Nhc, Npcs*Nwc*Nhc)
%       dAmatdw : (Npcs*Nwc*Nh, Npcs*Nwc*Nh) OR (Npcs*Nwc*Nhc, Npcs*Nwc*Nhc)
%       dAmatdxi: (Npcs*Nwc*Nh, Npcs*Nwc*Nh) OR (Npcs*Nwc*Nhc, Npcs*Nwc*Nhc)
%       Fv      : (Npcs*Nwc*Nh, 1) OR (Npcs*Nwc*Nhc, 1)
%       dFvdw   : (Npcs*Nwc*Nh, 1) OR (Npcs*Nwc*Nhc, 1)
%       dFvdxi  : (Npcs*Nwc*Nh, 1) OR (Npcs*Nwc*Nhc, 1)
%       JEV     : (array of structs) with fields
%           'Lj' , 'dLjdw' , 'dLjdxi'
%           'LjR', 'dLjRdw', 'dLjRdxi'
%           'Gj' , 'dGjdw' , 'dGjdxi'

    Npcs = length(pcs);
    Npts = pcs(end).irange(end);
    Nwc  = size(pcs(1).wcomps,1);  % Number of wave components
    Nh   = length(h);

    Amat = zeros(Npcs*Nwc*Nh);
    dAmatdw = zeros(Npcs*Nwc*Nh);
    dAmatdxi = zeros(Npcs*Nwc*Nh);
    Fv = zeros(Npcs*Nwc*Nh,1);
    dFvdw = zeros(Npcs*Nwc*Nh,1);
    dFvdxi = zeros(Npcs*Nwc*Nh,1);

    Rh = zeros(Npts*Nwc, Npcs*Nwc);
    dRhdw = zeros(Npts*Nwc, Npcs*Nwc);
    dRhdxi = zeros(Npts*Nwc, Npcs*Nwc);

    Ri = zeros(Npts*Nwc, 1);
    dRidw = zeros(Npts*Nwc, 1);
    dRidxi = zeros(Npts*Nwc, 1);

    % Setup Selector-Projectors for Nonlinearities
    nlis = find(arrayfun(@(j) ~isempty(j.nl), joints));
    nnl = length(nlis);
    Lj = cell(nnl,1);    dLjdw = cell(nnl,1);    dLjdxi = cell(nnl,1);
    LjR= cell(nnl,1);    dLjRdw= cell(nnl,1);    dLjRdxi= cell(nnl,1);
    Gj = cell(nnl,1);    dGjdw = cell(nnl,1);    dGjdxi = cell(nnl,1);
    for n=1:length(nlis)
        k = nlis(n);
        
        Lj{n} = zeros(joints(k).nld*Nh, Npcs*Nwc*Nh);
        dLjdw{n} = zeros(joints(k).nld*Nh, Npcs*Nwc*Nh);
        dLjdxi{n} = zeros(joints(k).nld*Nh, Npcs*Nwc*Nh);

        LjR{n} = zeros(joints(k).nld*Nh, 1);
        dLjRdw{n} = zeros(joints(k).nld*Nh, 1);
        dLjRdxi{n} = zeros(joints(k).nld*Nh, 1);

        Gj{n} = zeros(Npcs*Nwc*Nh, joints(k).nld*Nh);
        dGjdw{n} = zeros(Npcs*Nwc*Nh, joints(k).nld*Nh);
        dGjdxi{n} = zeros(Npcs*Nwc*Nh, joints(k).nld*Nh);
    end    
    
    % 1. Propagation Relationships
    for ih=1:Nh
        hstart_pt = (ih-1)*(Npts*Nwc);
        hstart_pc = (ih-1)*(Npcs*Nwc);

        for i=1:length(pcs)
            % Wavenumbers of the different components
            Ks = arrayfun(@(k) k.K(wxi(1)*h(ih), wxi(2)), Klib);
            dKdws = arrayfun(@(k) k.dKdw(wxi(1)*h(ih), wxi(2))*h(ih), Klib);
            dKdxis = arrayfun(@(k) k.dKdxi(wxi(1)*h(ih), wxi(2)), Klib);

            K = pcs(i).wcomps(:,1).*Ks(pcs(i).wcomps(:,2));
            dKdw = pcs(i).wcomps(:,1).*dKdws(pcs(i).wcomps(:,2));
            dKdxi = pcs(i).wcomps(:,1).*dKdxis(pcs(i).wcomps(:,2));

            % For first point
            si = hstart_pt*0 + ((pcs(i).irange(1)-1)+1-1)*Nwc+1;
            se = hstart_pt*0 + ((pcs(i).irange(1)-1)+1)*Nwc;

            qi = hstart_pc*0 + (i-1)*Nwc+1;
            qe = hstart_pc*0 + i*Nwc;

            inds = sub2ind([Npts Npcs]*Nwc, si:se, qi:qe);
            Rh(inds) = 1.0; % other terms are zero
            for n=2:pcs(i).N
                sim1 = si;
                sem1 = se;

                si = hstart_pt*0 + ((pcs(i).irange(1)-1)+n-1)*Nwc+1;
                se = hstart_pt*0 + ((pcs(i).irange(1)-1)+n)*Nwc;

                qi = hstart_pc*0 + (i-1)*Nwc+1;
                qe = hstart_pc*0 + i*Nwc;
                
                dx = abs((pcs(i).U(n)-pcs(i).U(n-1))*pcs(i).S);

                indsm1 = inds;
                inds = sub2ind([Npts Npcs]*Nwc, si:se, qi:qe);
                Rh(inds) = exp(K.'*dx).*Rh(indsm1);
                dRhdw(inds) = Rh(inds).*dKdw.'*dx + ...
                    exp(K.'*dx).*dRhdw(indsm1);
                dRhdxi(inds) = Rh(inds).*dKdxi.'*dx + ...
                    exp(K.'*dx).*dRhdxi(indsm1);

                Ri(si:se) = exp(K*dx).*Ri(sim1:sem1);
                dRidw(si:se) = Ri(si:se).*dKdw*dx + ...
                    exp(K*dx).*dRidw(sim1:sem1);
                dRidxi(si:se) = Ri(si:se).*dKdxi*dx + ...
                    exp(K*dx).*dRidxi(sim1:sem1);

                % Check if this is an excitation point
                ei = find(pcs(i).exci==n-1);
                if ~isempty(ei)
                    if length(ei)>1
                        error('More than 1 excitation location detected. Check pcs.');
                    end

                    if pcs(i).excnh(ei)==h(ih)  % SHOULD THESE BE NEGATIVE?
                        Ri(si:se) = Ri(si:se) + pcs(i).exccofs{ei}{1}(h(ih)*wxi(1), wxi(2));
                        dRidw(si:se) = dRidw(si:se) + pcs(i).exccofs{ei}{2}(h(ih)*wxi(1), wxi(2))*h(ih);
                        dRidxi(si:se) = dRidxi(si:se) + pcs(i).exccofs{ei}{2}(h(ih)*wxi(1), wxi(2));
                    end
                end
            end
        end
        
        hinds = hstart_pc + (1:Npcs*Nwc);
        % 2. Boundary Conditions
        for n=1:length(bcs)
            inds = hstart_pt*0 + (bcs(n).i-1)*Nwc+(1:Nwc);
            
            % RESTRICT THE COLUMNS OF THESE EQNS TO THE CURRENT HARMONIC ONLY.
            Amat(hstart_pc+n, hinds) = bcs(n).cofs(h(ih)*wxi(1), wxi(2))*Rh(inds,hinds-hstart_pc);
            dAmatdw(hstart_pc+n, hinds) = bcs(n).dcofsdw(h(ih)*wxi(1), wxi(2))*Rh(inds,hinds-hstart_pc)*h(ih) + ...
                bcs(n).cofs(h(ih)*wxi(1), wxi(2))*dRhdw(inds,hinds-hstart_pc);
            dAmatdxi(hstart_pc+n, hinds) = bcs(n).dcofsdxi(h(ih)*wxi(1), wxi(2))*Rh(inds,hinds-hstart_pc) + ...
                bcs(n).cofs(h(ih)*wxi(1), wxi(2))*dRhdxi(inds,hinds-hstart_pc);

            Fv(hstart_pc+n) = bcs(n).cofs(h(ih)*wxi(1), wxi(2))*Ri(inds);
            dFvdw(hstart_pc+n) = bcs(n).dcofsdw(h(ih)*wxi(1), wxi(2))*Ri(inds)*h(ih) + ...
                bcs(n).cofs(h(ih)*wxi(1), wxi(2))*dRidw(inds);
            dFvdxi(hstart_pc+n) = bcs(n).dcofsdxi(h(ih)*wxi(1), wxi(2))*Ri(inds) + ...
                bcs(n).cofs(h(ih)*wxi(1), wxi(2))*dRidxi(inds);

            % RHS
            if ~isempty(bcs(n).rih) && bcs(n).rih==h(ih)
                Fv(hstart_pc+n) = Fv(hstart_pc+n) + bcs(n).rhs(h(ih)*wxi(1), wxi(2));
                dFvdw(hstart_pc+n) = dFvdw(hstart_pc+n) + bcs(n).drhsdw(h(ih)*wxi(1), wxi(2))*h(ih);
                dFvdxi(hstart_pc+n) = dFvdxi(hstart_pc+n) + bcs(n).drhsdxi(h(ih)*wxi(1), wxi(2));
            end
        end

        % 3. Joints
        ktn = hstart_pc + length(bcs);
        for n=1:length(joints)
            switch joints(n).type
                case {1, 2}  % Binary Connection
                    inds = hstart_pt*0 + [(joints(n).i-1)*Nwc+(1:Nwc) (joints(n).j-1)*Nwc+(1:Nwc)];
                    
                    % RESTRICT THE COLUMNS OF THESE EQNS TO THE CURRENT HARMONIC ONLY.
                    Amat(ktn+(n-1)*Nwc+(1:Nwc), hinds) = joints(n).cofs(h(ih)*wxi(1), wxi(2))*Rh(inds,hinds-hstart_pc);
                    dAmatdw(ktn+(n-1)*Nwc+(1:Nwc), hinds) = joints(n).dcofsdw(h(ih)*wxi(1), wxi(2))*Rh(inds,hinds-hstart_pc)*h(ih) + ...
                        joints(n).cofs(h(ih)*wxi(1), wxi(2))*dRhdw(inds,hinds-hstart_pc);
                    dAmatdxi(ktn+(n-1)*Nwc+(1:Nwc), hinds) = joints(n).dcofsdxi(h(ih)*wxi(1), wxi(2))*Rh(inds,hinds-hstart_pc) + ...
                        joints(n).cofs(h(ih)*wxi(1), wxi(2))*dRhdxi(inds,hinds-hstart_pc);

                    Fv(ktn+(n-1)*Nwc+(1:Nwc)) = joints(n).cofs(h(ih)*wxi(1), wxi(2))*Ri(inds);
                    dFvdw(ktn+(n-1)*Nwc+(1:Nwc)) = joints(n).dcofsdw(h(ih)*wxi(1), wxi(2))*Ri(inds)*h(ih) + ...
                        joints(n).cofs(h(ih)*wxi(1), wxi(2))*dRidw(inds);
                    dFvdxi(ktn+(n-1)*Nwc+(1:Nwc)) = joints(n).dcofsdxi(h(ih)*wxi(1), wxi(2))*Ri(inds) + ...
                        joints(n).cofs(h(ih)*wxi(1), wxi(2))*dRidxi(inds);
                otherwise
                    error('Needs to be implemented still.');
            end
        end
        % 3a. Nonlinear Joints
        for n=1:nnl
            k = nlis(n);
            nld = joints(k).nld;

            switch joints(k).type
                case 2  
                    % Choosing NL displacement
                    inds = hstart_pt*0 + [(joints(k).i-1)*Nwc+(1:Nwc) (joints(k).j-1)*Nwc+(1:Nwc)];

                    Lj{n}((ih-1)*nld+(1:nld), hinds) = joints(k).nldcofs(h(ih)*wxi(1), wxi(2))*Rh(inds,hinds-hstart_pc);
                    dLjdw{n}((ih-1)*nld+(1:nld), hinds) = joints(k).dnldcofsdw(h(ih)*wxi(1), wxi(2))*Rh(inds,hinds-hstart_pc)*h(ih) + ...
                        joints(k).nldcofs(h(ih)*wxi(1), wxi(2))*dRhdw(inds,hinds-hstart_pc);
                    dLjdxi{n}((ih-1)*nld+(1:nld), hinds) = joints(k).dnldcofsdxi(h(ih)*wxi(1), wxi(2))*Rh(inds,hinds-hstart_pc) + ...
                        joints(k).nldcofs(h(ih)*wxi(1), wxi(2))*dRhdxi(inds,hinds-hstart_pc);

                    LjR{n}((ih-1)*nld+(1:nld)) = joints(k).nldcofs(h(ih)*wxi(1), wxi(2))*Ri(inds);
                    dLjRdw{n}((ih-1)*nld+(1:nld)) = joints(k).dnldcofsdw(h(ih)*wxi(1), wxi(2))*Ri(inds)*h(ih) + ...
                        joints(k).nldcofs(h(ih)*wxi(1), wxi(2))*dRidw(inds);
                    dLjRdxi{n}((ih-1)*nld+(1:nld)) = joints(k).dnldcofsdxi(h(ih)*wxi(1), wxi(2))*Ri(inds) + ...
                        joints(k).nldcofs(h(ih)*wxi(1), wxi(2))*dRidxi(inds);

                    % Putting NL Force
                    kinds = ktn+(k-1)*Nwc+(1:Nwc);

                    Gj{n}(kinds, (ih-1)*nld+(1:nld)) = joints(k).nlfcofs(h(ih)*wxi(1), wxi(2));
                    dGjdw{n}(kinds, (ih-1)*nld+(1:nld)) = joints(k).dnlfcofsdw(h(ih)*wxi(1), wxi(2))*h(ih);
                    dGjdxi{n}(kinds, (ih-1)*nld+(1:nld)) = joints(k).dnlfcofsdxi(h(ih)*wxi(1), wxi(2));
                otherwise
                    error('Needs to be implemented still.');
            end
            if ~isempty(joints(k).rih) && joints(k).rih==h(ih)
                Fv(kinds) = Fv(kinds) + joints(k).rhs(h(ih)*wxi(1), wxi(2));
                dFvdw(kinds) = dFvdw(kinds) + joints(k).drhsdw(h(ih)*wxi(1), wxi(2))*h(ih);
                dFvdxi(kinds) = dFvdxi(kinds) + joints(k).drhsdxi(h(ih)*wxi(1), wxi(2));
            end
        end
    end
    JEV = struct('Lj', Lj, 'dLjdw', dLjdw, 'dLjdxi', dLjdxi, ...
        'LjR', LjR, 'dLjRdw', dLjRdw, 'dLjRdxi', dLjRdxi, ...
        'Gj', Gj, 'dGjdw', dGjdw, 'dGjdxi', dGjdxi);

    %% Convert to Fully Real Representation
    if length(varargin)>=1 && varargin{1}=='r'
        Nhc = sum((h==0)+2*(h~=0));
        [zinds,hinds,rinds0,rinds,iinds] = HINDS(Npcs*Nwc, h);

        Amatc = Amat;
        dAmatdwc = dAmatdw;
        dAmatdxic = dAmatdxi;
    
        Amat = zeros(Npcs*Nwc*Nhc);
        dAmatdw = zeros(Npcs*Nwc*Nhc);
        dAmatdxi = zeros(Npcs*Nwc*Nhc);
    
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

        Fv = zeros(Npcs*Nwc*Nhc,1);
        dFvdw = zeros(Npcs*Nwc*Nhc,1);
        dFvdxi = zeros(Npcs*Nwc*Nhc,1);
        Fv([rinds0 rinds iinds]) = [real(Fvc(zinds)); real(Fvc(hinds)); imag(Fvc(hinds))];
        dFvdw([rinds0 rinds iinds]) = [real(dFvdwc(zinds)); real(dFvdwc(hinds)); imag(dFvdwc(hinds))];
        dFvdxi([rinds0 rinds iinds]) = [real(dFvdxic(zinds)); real(dFvdxic(hinds)); imag(dFvdxic(hinds))];

        % Nonlinear Selector-Projectors
        JEVc = JEV;
        JEV = struct('Lj', cell(nnl,1), 'dLjdw', cell(nnl,1), 'dLjdxi', cell(nnl,1), ...
            'LjR', cell(nnl,1), 'dLjRdw', cell(nnl,1), 'dLjRdxi', cell(nnl,1), ...
            'Gj', cell(nnl,1), 'dGjdw', cell(nnl,1), 'dGjdxi', cell(nnl,1));

        for n=1:nnl
            k = nlis(n);
            nld = joints(k).nld;

            [nlzinds,nlhinds,nlrinds0,nlrinds,nliinds] = HINDS(nld, h);
            JEV(n).Lj = zeros(nld*Nhc, Npcs*Nwc*Nhc);
            JEV(n).dLjdw = zeros(nld*Nhc, Npcs*Nwc*Nhc);
            JEV(n).dLjdxi = zeros(nld*Nhc, Npcs*Nwc*Nhc);

            JEV(n).LjR = zeros(nld*Nhc, 1);
            JEV(n).dLjRdw = zeros(nld*Nhc, 1);
            JEV(n).dLjRdxi = zeros(nld*Nhc, 1);
            
            JEV(n).Gj = zeros(Npcs*Nwc*Nhc, nld*Nhc);
            JEV(n).dGjdw = zeros(Npcs*Nwc*Nhc, nld*Nhc);
            JEV(n).dGjdxi = zeros(Npcs*Nwc*Nhc, nld*Nhc);

            tmp = [JEVc(n).Lj(nlhinds, hinds) 1j*JEVc(n).Lj(nlhinds, hinds)];
            JEV(n).Lj(nlrinds0, rinds0) = JEVc(n).Lj(nlzinds, zinds);
            JEV(n).Lj([nlrinds nliinds], [rinds iinds]) = [real(tmp);imag(tmp)];
            tmp = [JEVc(n).dLjdw(nlhinds, hinds) 1j*JEVc(n).dLjdw(nlhinds, hinds)];
            JEV(n).dLjdw(nlrinds0, rinds0) = JEVc(n).dLjdw(nlzinds, zinds);
            JEV(n).dLjdw([nlrinds nliinds], [rinds iinds]) = [real(tmp);imag(tmp)];
            tmp = [JEVc(n).dLjdxi(nlhinds, hinds) 1j*JEVc(n).dLjdxi(nlhinds, hinds)];
            JEV(n).dLjdxi(nlrinds0, rinds0) = JEVc(n).dLjdxi(nlzinds, zinds);
            JEV(n).dLjdxi([nlrinds nliinds], [rinds iinds]) = [real(tmp);imag(tmp)];

            JEV(n).LjR([nlrinds0 nlrinds nliinds]) = [JEVc(n).LjR(nlzinds); real(JEVc(n).LjR(nlhinds)); imag(JEVc(n).LjR(nlhinds))];
            JEV(n).dLjRdw([nlrinds0 nlrinds nliinds]) = [JEVc(n).dLjRdw(nlzinds); real(JEVc(n).dLjRdw(nlhinds)); imag(JEVc(n).dLjRdw(nlhinds))];
            JEV(n).dLjRdxi([nlrinds0 nlrinds nliinds]) = [JEVc(n).dLjRdxi(nlzinds); real(JEVc(n).dLjRdxi(nlhinds)); imag(JEVc(n).dLjRdxi(nlhinds))];

            tmp = [JEVc(n).Gj(hinds,nlhinds) 1j*JEVc(n).Gj(hinds,nlhinds)];
            JEV(n).Gj(rinds0,nlrinds0) = JEVc(n).Gj(zinds,nlzinds);
            JEV(n).Gj([rinds iinds],[nlrinds nliinds]) = [real(tmp);imag(tmp)];
            tmp = [JEVc(n).dGjdw(hinds,nlhinds) 1j*JEVc(n).dGjdw(hinds,nlhinds)];
            JEV(n).dGjdw(rinds0,nlrinds0) = JEVc(n).dGjdw(zinds,nlzinds);
            JEV(n).dGjdw([rinds iinds],[nlrinds nliinds]) = [real(tmp);imag(tmp)];
            tmp = [JEVc(n).dGjdxi(hinds,nlhinds) 1j*JEVc(n).dGjdxi(hinds,nlhinds)];
            JEV(n).dGjdxi(rinds0,nlrinds0) = JEVc(n).dGjdxi(zinds,nlzinds);
            JEV(n).dGjdxi([rinds iinds],[nlrinds nliinds]) = [real(tmp);imag(tmp)];
        end
    end    
end