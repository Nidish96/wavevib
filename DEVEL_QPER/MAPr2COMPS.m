function [Rh, dRhdw, dRhdxi, Ri, dRidw, dRidxi] = MAPr2COMPS(wxi, h, pcs, Klib)
%MAP22COMPS returns mapping matrices to map the reduced parameterization
%(one set of wave components per piece) to the full set (one set at each
%provided node).
%       Mapping is implemented such that
%           Afull = Rh*Ared - Ri 
%
%   USAGE:
%       [Rh, dRhdw, dRhdxi, Ri, dRidw, dRidxi] = MAPr2COMPS(wxi, h, pcs, Klib);
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
%       Klib    : (array of structs) function handles of dispersion
%           relationships in the form of K=f(w,xi) with fields
%           'K', 'dKdw', 'dKdxi'
%   OUTPUTS:
%       Rh      : (Npts*Nwc*Nh, Npcs*Nwc*Nh) Complex matrix
%       dRhdw   : (Npts*Nwc*Nh, Npcs*Nwc*Nh, Nc)
%       dRhdxi  : (Npts*Nwc*Nh, Npcs*Nwc*Nh, Npar)
%       Ri      : (Npts*Nwc*Nh, 1) Complex vector
%       dRidw   : (Npts*Nwc*Nh, Nc)
%       dRidxi  : (Npts*Nwc*Nh, Npar)

    Npcs = length(pcs);
    Npts = pcs(end).irange(end);
    Nwc  = size(pcs(1).wcomps,1);  % Number of wave components
    Nh   = size(h,1);  % Number of harmonics
    Nc   = size(h,2);  % Number of incommensurate frequency components
    Npar = length(wxi)-Nc;

    Rh = zeros(Npts*Nwc*Nh, Npcs*Nwc*Nh);
    dRhdw = zeros((Npts*Nwc*Nh)*(Npcs*Nwc*Nh), Nc);
    dRhdxi = zeros((Npts*Nwc*Nh)*(Npcs*Nwc*Nh), Npar);

    Ri = zeros(Npts*Nwc*Nh, 1);
    dRidw = zeros(Npts*Nwc*Nh, Nc);
    dRidxi = zeros(Npts*Nwc*Nh, Npar);

    for ih=1:Nh
        hstart_pt = (ih-1)*(Npts*Nwc);
        hstart_pc = (ih-1)*(Npcs*Nwc);
        
        for i=1:length(pcs)
            % Wavenumbers of the different components
            Ks = arrayfun(@(k) k.K(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))), Klib);  %(nk,1)
            dKdws = cell2mat(arrayfun(@(k) k.dKdw(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)))*h(ih, :), Klib, 'UniformOutput', false));  %(nk,Nc)
            dKdxis = arrayfun(@(k) k.dKdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))), Klib);  %(nk,Npar)

            K = pcs(i).wcomps(:,1).*Ks(pcs(i).wcomps(:,2));  %(Nwc,1)
            dKdw = pcs(i).wcomps(:,1).*dKdws(pcs(i).wcomps(:,2),:);  %(Nwc,Nc)
            dKdxi = pcs(i).wcomps(:,1).*dKdxis(pcs(i).wcomps(:,2),:);  %(Nwc,Npar)

            % For first point
            si = hstart_pt + ((pcs(i).irange(1)-1)+1-1)*Nwc+1;
            se = hstart_pt + ((pcs(i).irange(1)-1)+1)*Nwc;

            qi = hstart_pc + (i-1)*Nwc+1;
            qe = hstart_pc + i*Nwc;

            inds = sub2ind([Npts Npcs]*Nwc*Nh, (si:se), (qi:qe));  %% Fix indices
            Rh(inds) = 1.0; % other terms are zero            
            Ri(si:se) = 0.0;
            for n=2:pcs(i).N
                sim1 = si;
                sem1 = se;

                si = hstart_pt + ((pcs(i).irange(1)-1)+n-1)*Nwc+1;
                se = hstart_pt + ((pcs(i).irange(1)-1)+n)*Nwc;

                qi = hstart_pc + (i-1)*Nwc+1;
                qe = hstart_pc + i*Nwc;
                
                dx = abs((pcs(i).U(n)-pcs(i).U(n-1))*pcs(i).S);

                indsm1 = inds;
                inds = sub2ind([Npts Npcs]*Nwc*Nh, si:se, qi:qe);
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
    end
    dRhdw = reshape(dRhdw, Npts*Nwc*Nh, Npcs*Nwc*Nh, Nc);
    dRhdxi = reshape(dRhdxi, Npts*Nwc*Nh, Npcs*Nwc*Nh, Npar);
end