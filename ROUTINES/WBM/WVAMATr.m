function [Amat, dAmatdw, dAmatdxi, Fv, dFvdw, dFvdxi, JEV] = WVAMATr(wxi, h, pcs, bcs, joints, Klib, varargin)
%WVAMATr returns the linear "A" matrix and "F" vector for the wave-based
%model using reduced parameterization. Also returned are the Jacobians w.r.t.
%w (frequency) and xi (parameter) and selector-projector matrices for
%evaluating the nonlinearities.
% Returns in either complex or real-imaginary representation.
%   This routine is for the general quasi-periodic case (with multiple
%   incommensurate frequency components).
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
%       wxi     : (Nc+Npar,1) vector of (frequency,parameter)
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
    Nh   = size(h,1);  % Number of harmonics
    Nc   = size(h,2);  % Number of incommensurate frequency components
    Npar = length(wxi)-Nc;

    Amat = zeros(Npcs*Nwc*Nh);
    dAmatdw = zeros(Npcs*Nwc*Nh,Npcs*Nwc*Nh,Nc);
    dAmatdxi = zeros(Npcs*Nwc*Nh,Npcs*Nwc*Nh,Npar);
    Fv = zeros(Npcs*Nwc*Nh,1);
    dFvdw = zeros(Npcs*Nwc*Nh,Nc);
    dFvdxi = zeros(Npcs*Nwc*Nh,Npar);

    Rh = zeros(Npts*Nwc, Npcs*Nwc);
    dRhdw = zeros((Npts*Nwc)*(Npcs*Nwc), Nc);
    dRhdxi = zeros((Npts*Nwc)*(Npcs*Nwc), Npar);

    Ri = zeros(Npts*Nwc, 1);
    dRidw = zeros(Npts*Nwc, Nc);
    dRidxi = zeros(Npts*Nwc, Npar);

    % Setup Selector-Projectors for Nonlinearities
    nlis = find(arrayfun(@(j) ~isempty(j.nl), joints));
    nnl = length(nlis);
    Lj = cell(nnl,1);    dLjdw = cell(nnl,1);    dLjdxi = cell(nnl,1);
    LjR= cell(nnl,1);    dLjRdw= cell(nnl,1);    dLjRdxi= cell(nnl,1);
    Gj = cell(nnl,1);    dGjdw = cell(nnl,1);    dGjdxi = cell(nnl,1);
    for n=1:length(nlis)
        k = nlis(n);
        
        Lj{n} = zeros(joints(k).nld*Nh, Npcs*Nwc*Nh);
        dLjdw{n} = zeros(joints(k).nld*Nh, Npcs*Nwc*Nh, Nc);
        dLjdxi{n} = zeros(joints(k).nld*Nh, Npcs*Nwc*Nh, Npar);

        LjR{n} = zeros(joints(k).nld*Nh, 1);
        dLjRdw{n} = zeros(joints(k).nld*Nh, Nc);
        dLjRdxi{n} = zeros(joints(k).nld*Nh, Npar);

        Gj{n} = zeros(Npcs*Nwc*Nh, joints(k).nld*Nh);
        dGjdw{n} = zeros(Npcs*Nwc*Nh, joints(k).nld*Nh, Nc);
        dGjdxi{n} = zeros(Npcs*Nwc*Nh, joints(k).nld*Nh, Npar);
    end    
    
    % 1. Propagation Relationships
    for ih=1:Nh
        hstart_pt = (ih-1)*(Npts*Nwc);
        hstart_pc = (ih-1)*(Npcs*Nwc);

        Rh = reshape(Rh, Npts*Nwc, Npcs*Nwc);
        dRhdw = reshape(dRhdw, (Npts*Nwc)*(Npcs*Nwc), Nc);
        dRhdxi = reshape(dRhdxi, (Npts*Nwc)*(Npcs*Nwc), Npar);

        for i=1:length(pcs)
            % Wavenumbers of the different components
            Ks = arrayfun(@(k) k.K(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))), Klib);  %(nk,1)
            dKdws = cell2mat(arrayfun(@(k) k.dKdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*h(ih, :), Klib, 'UniformOutput', false));  %(nk,Nc)
            dKdxis = arrayfun(@(k) k.dKdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))), Klib);  %(nk,Npar)

            K = pcs(i).wcomps(:,1).*Ks(pcs(i).wcomps(:,2));  %(Nwc,1)
            dKdw = pcs(i).wcomps(:,1).*dKdws(pcs(i).wcomps(:,2),:);  %(Nwc,Nc)
            dKdxi = pcs(i).wcomps(:,1).*dKdxis(pcs(i).wcomps(:,2),:);  %(Nwc,Npar)

            % For first point
            si = hstart_pt*0 + ((pcs(i).irange(1)-1)+1-1)*Nwc+1;
            se = hstart_pt*0 + ((pcs(i).irange(1)-1)+1)*Nwc;

            qi = hstart_pc*0 + (i-1)*Nwc+1;
            qe = hstart_pc*0 + i*Nwc;

            inds = sub2ind([Npts Npcs]*Nwc, si:se, qi:qe);  %% Fix indices
            Rh(inds) = 1.0; % other terms are zero
            Ri(si:se) = 0.0;
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
                Rh(inds(:)) = exp(K*dx).*Rh(indsm1(:));
                dRhdw(inds,:) = Rh(inds(:)).*dKdw*dx + ...
                    exp(K*dx).*dRhdw(indsm1,:);
                dRhdxi(inds,:) = Rh(inds(:)).*dKdxi*dx + ...
                    exp(K*dx).*dRhdxi(indsm1,:);

                Ri(si:se) = exp(K*dx).*Ri(sim1:sem1);
                dRidw(si:se,:) = Ri(si:se).*dKdw*dx + ...
                    exp(K*dx).*dRidw(sim1:sem1,:);
                dRidxi(si:se,:) = Ri(si:se).*dKdxi*dx + ...
                    exp(K*dx).*dRidxi(sim1:sem1,:);

                % Check if this is an excitation point
                ei = find(pcs(i).exci==n-1);
                if ~isempty(ei)
                    if length(ei)>1
                        error('More than 1 excitation location detected. Check pcs.');
                    end

                    if all(pcs(i).excnh(ei,:)==h(ih,:))
                        Ri(si:se) = Ri(si:se) + pcs(i).exccofs{ei}{1}(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));
                        dRidw(si:se,:) = dRidw(si:se,:) + pcs(i).exccofs{ei}{2}(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))).*h(ih,:);
                        dRidxi(si:se,:) = dRidxi(si:se,:) + pcs(i).exccofs{ei}{2}(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));
                    end
                end
            end
        end
        dRhdw = reshape(dRhdw, Npts*Nwc, Npcs*Nwc, Nc);
        dRhdxi = reshape(dRhdxi, Npts*Nwc, Npcs*Nwc, Npar);
        
        hinds = hstart_pc + (1:Npcs*Nwc);
        % 2. Boundary Conditions
        nofs = 0;
        for n=1:length(bcs)
            nofs = [0 bcs.nof];

            inds = hstart_pt*0 + (bcs(n).i-1)*Nwc+(1:Nwc);
            
            % RESTRICT THE COLUMNS OF THESE EQNS TO THE CURRENT HARMONIC ONLY.
            Amat(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1)), hinds) = bcs(n).cofs(h(ih, :)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Rh(inds,hinds-hstart_pc);
            dAmatdw(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1)), hinds, :) = bcs(n).dcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Rh(inds,hinds-hstart_pc).*permute(h(ih,:),[1 3 2]) + ...
                reshape(cell2mat(arrayfun(@(a) bcs(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRhdw(inds,hinds-hstart_pc,a), 1:Nc, 'UniformOutput', false)), nofs(n+1), Npcs*Nwc, Nc);
            dAmatdxi(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1)), hinds, :) = bcs(n).dcofsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Rh(inds,hinds-hstart_pc) + ...  %%%%
                reshape(cell2mat(arrayfun(@(a) bcs(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRhdxi(inds,hinds-hstart_pc,a), 1:Npar, 'UniformOutput', false)), nofs(n+1), Npcs*Nwc, Npar);

            Fv(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1))) = bcs(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Ri(inds);
            dFvdw(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1)),:) = bcs(n).dcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Ri(inds).*h(ih,:) + ...
                bcs(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRidw(inds,:);
            dFvdxi(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1)),:) = bcs(n).dcofsdxi(h(ih, :)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Ri(inds) + ...
                bcs(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRidxi(inds, :);

            % RHS
            if ~isempty(bcs(n).rih) && all(bcs(n).rih==h(ih,:))
                Fv(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1))) = Fv(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1))) + bcs(n).rhs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));
                dFvdw(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1)),:) = dFvdw(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1)),:) + bcs(n).drhsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))).*h(ih,:);
                dFvdxi(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1)),:) = dFvdxi(hstart_pc+sum(nofs(1:n))+(1:nofs(n+1)),:) + bcs(n).drhsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));
            end
        end

        % 3. Joints
        ktn = hstart_pc + sum(nofs);
        for n=1:length(joints)
            nofs = [0 joints.nof];
            switch joints(n).type
                case {1, 2}  % Binary Connection
                    inds = hstart_pt*0 + [(joints(n).i-1)*Nwc+(1:Nwc) (joints(n).j-1)*Nwc+(1:Nwc)];
                otherwise
                    inds = hstart_pt*0 + cell2mat(arrayfun(@(i) (i-1)*Nwc+(1:Nwc), joints(n).is, 'UniformOutput', false));
            end

            % RESTRICT THE COLUMNS OF THESE EQNS TO THE CURRENT HARMONIC ONLY.
            Amat(ktn+sum(nofs(1:n))+(1:nofs(n+1)), hinds) = joints(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Rh(inds,hinds-hstart_pc);
            dAmatdw(ktn+sum(nofs(1:n))+(1:nofs(n+1)), hinds, :) = joints(n).dcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Rh(inds,hinds-hstart_pc).*permute(h(ih,:),[1 3 2]) + ...
                reshape(cell2mat(arrayfun(@(a) joints(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRhdw(inds,hinds-hstart_pc,a), 1:Nc, 'UniformOutput',false)), nofs(n+1), Npcs*Nwc, Nc);
            dAmatdxi(ktn+sum(nofs(1:n))+(1:nofs(n+1)), hinds, :) = joints(n).dcofsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Rh(inds,hinds-hstart_pc) + ...  %%%%
                reshape(cell2mat(arrayfun(@(a) joints(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRhdxi(inds,hinds-hstart_pc,a), 1:Npar, 'UniformOutput', false)), nofs(n+1), Npcs*Nwc, Npar);

            Fv(ktn+sum(nofs(1:n))+(1:nofs(n+1))) = joints(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Ri(inds);
            dFvdw(ktn+sum(nofs(1:n))+(1:nofs(n+1)),:) = joints(n).dcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Ri(inds).*h(ih,:) + ...
                joints(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRidw(inds,:);
            dFvdxi(ktn+sum(nofs(1:n))+(1:nofs(n+1))) = joints(n).dcofsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Ri(inds) + ...
                joints(n).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRidxi(inds,:);
        end
        % 3a. Nonlinear Joints
        for n=1:nnl
            k = nlis(n);
            nld = joints(k).nld;

            switch joints(k).type
                case 2  
                    % Choosing NL displacement
                    inds = hstart_pt*0 + [(joints(k).i-1)*Nwc+(1:Nwc) (joints(k).j-1)*Nwc+(1:Nwc)];
                otherwise
                    inds = hstart_pt*0 + cell2mat(arrayfun(@(i) (i-1)*Nwc+(1:Nwc), joints(k).is, 'UniformOutput', false));
            end
            Lj{n}((ih-1)*nld+(1:nld), hinds) = joints(k).nldcofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Rh(inds,hinds-hstart_pc);
            dLjdw{n}((ih-1)*nld+(1:nld), hinds, :) = joints(k).dnldcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Rh(inds,hinds-hstart_pc).*permute(h(ih,:), [1 3 2]) + ...
                reshape(cell2mat(arrayfun(@(a) joints(k).nldcofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRhdw(inds,hinds-hstart_pc,a), 1:Nc, 'UniformOutput', false)), nld, Npcs*Nwc, Nc);
            dLjdxi{n}((ih-1)*nld+(1:nld), hinds, :) = joints(k).dnldcofsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Rh(inds,hinds-hstart_pc) + ...  %%%%
                reshape(cell2mat(arrayfun(@(a) joints(k).nldcofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRhdxi(inds,hinds-hstart_pc,a), 1:Npar, 'UniformOutput', false)), nld, Npcs*Nwc, Npar);

            LjR{n}((ih-1)*nld+(1:nld)) = joints(k).nldcofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Ri(inds);
            dLjRdw{n}((ih-1)*nld+(1:nld),:) = joints(k).dnldcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Ri(inds).*h(ih,:) + ...
                joints(k).nldcofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRidw(inds,:);
            dLjRdxi{n}((ih-1)*nld+(1:nld),:) = joints(k).dnldcofsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*Ri(inds) + ...
                joints(k).nldcofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*dRidxi(inds,:);

            % Putting NL Force
            kinds = ktn+sum(nofs(1:k))+(1:nofs(k+1));
%             kinds = ktn+(k-1)*Nwc+(1:Nwc);            

            Gj{n}(kinds, (ih-1)*nld+(1:nld)) = joints(k).nlfcofs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));
            dGjdw{n}(kinds, (ih-1)*nld+(1:nld), :) = joints(k).dnlfcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))).*permute(h(ih,:), [1 3 2]);
            dGjdxi{n}(kinds, (ih-1)*nld+(1:nld), :) = joints(k).dnlfcofsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));            

            if ~isempty(joints(k).rih) && all(joints(k).rih==h(ih,:))
                Fv(kinds) = Fv(kinds) + joints(k).rhs(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));
                dFvdw(kinds,:) = dFvdw(kinds,:) + joints(k).drhsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))).*h(ih,:);
                dFvdxi(kinds,:) = dFvdxi(kinds,:) + joints(k).drhsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));
            end
        end
    end
    JEV = struct('Lj', Lj, 'dLjdw', dLjdw, 'dLjdxi', dLjdxi, ...
        'LjR', LjR, 'dLjRdw', dLjRdw, 'dLjRdxi', dLjRdxi, ...
        'Gj', Gj, 'dGjdw', dGjdw, 'dGjdxi', dGjdxi);

    %% Convert to Fully Real Representation
    if length(varargin)>=1 && varargin{1}=='r'
        Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
        [zinds,hinds,rinds0,rinds,iinds] = HINDS(Npcs*Nwc, h);

        Amatc = Amat;
        dAmatdwc = dAmatdw;
        dAmatdxic = dAmatdxi;
    
        Amat = zeros(Npcs*Nwc*Nhc);
        dAmatdw = zeros(Npcs*Nwc*Nhc,Npcs*Nwc*Nhc,Nc);
        dAmatdxi = zeros(Npcs*Nwc*Nhc,Npcs*Nwc*Nhc,Npar);
    
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

        Fv = zeros(Npcs*Nwc*Nhc,1);
        dFvdw = zeros(Npcs*Nwc*Nhc,Nc);
        dFvdxi = zeros(Npcs*Nwc*Nhc,Npar);
        Fv([rinds0 rinds iinds]) = [real(Fvc(zinds)); real(Fvc(hinds)); imag(Fvc(hinds))];
        dFvdw([rinds0 rinds iinds], :) = [real(dFvdwc(zinds, :)); real(dFvdwc(hinds, :)); imag(dFvdwc(hinds, :))];
        dFvdxi([rinds0 rinds iinds], :) = [real(dFvdxic(zinds, :)); real(dFvdxic(hinds, :)); imag(dFvdxic(hinds, :))];

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
            JEV(n).dLjdw = zeros(nld*Nhc, Npcs*Nwc*Nhc, Nc);
            JEV(n).dLjdxi = zeros(nld*Nhc, Npcs*Nwc*Nhc, Npar);

            JEV(n).LjR = zeros(nld*Nhc, 1);
            JEV(n).dLjRdw = zeros(nld*Nhc, Nc);
            JEV(n).dLjRdxi = zeros(nld*Nhc, Npar);
            
            JEV(n).Gj = zeros(Npcs*Nwc*Nhc, nld*Nhc);
            JEV(n).dGjdw = zeros(Npcs*Nwc*Nhc, nld*Nhc, Nc);
            JEV(n).dGjdxi = zeros(Npcs*Nwc*Nhc, nld*Nhc, Npar);

            tmp = [JEVc(n).Lj(nlhinds, hinds) 1j*JEVc(n).Lj(nlhinds, hinds)];
            JEV(n).Lj(nlrinds0, rinds0) = JEVc(n).Lj(nlzinds, zinds);
            JEV(n).Lj([nlrinds nliinds], [rinds iinds]) = [real(tmp);imag(tmp)];
            tmp = [JEVc(n).dLjdw(nlhinds, hinds, :) 1j*JEVc(n).dLjdw(nlhinds, hinds, :)];
            JEV(n).dLjdw(nlrinds0, rinds0, :) = JEVc(n).dLjdw(nlzinds, zinds, :);
            JEV(n).dLjdw([nlrinds nliinds], [rinds iinds], :) = [real(tmp);imag(tmp)];
            tmp = [JEVc(n).dLjdxi(nlhinds, hinds, :) 1j*JEVc(n).dLjdxi(nlhinds, hinds, :)];
            JEV(n).dLjdxi(nlrinds0, rinds0, :) = JEVc(n).dLjdxi(nlzinds, zinds, :);
            JEV(n).dLjdxi([nlrinds nliinds], [rinds iinds], :) = [real(tmp);imag(tmp)];

            JEV(n).LjR([nlrinds0 nlrinds nliinds]) = [JEVc(n).LjR(nlzinds); real(JEVc(n).LjR(nlhinds)); imag(JEVc(n).LjR(nlhinds))];
            JEV(n).dLjRdw([nlrinds0 nlrinds nliinds], :) = [JEVc(n).dLjRdw(nlzinds, :); real(JEVc(n).dLjRdw(nlhinds, :)); imag(JEVc(n).dLjRdw(nlhinds, :))];
            JEV(n).dLjRdxi([nlrinds0 nlrinds nliinds], :) = [JEVc(n).dLjRdxi(nlzinds, :); real(JEVc(n).dLjRdxi(nlhinds, :)); imag(JEVc(n).dLjRdxi(nlhinds, :))];

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
