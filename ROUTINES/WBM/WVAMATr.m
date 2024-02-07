function [Amatr, dAmatdwr, dAmatdxir, Fvr, dFvdwr, dFvdxir, JEVr, RECOV] = WVAMATr(wxi, h, pcs, bcs, joints, Klib, varargin)

    Npts = pcs(end).irange(end);
    Nc = size(h,2);
    Npar = length(wxi)-Nc;
    Nwc = size(pcs(1).wcomps, 1);

    Nh = length(h);
    
    % Compute full Amat
    [Amat, dAmatdw, dAmatdxi, Fv, dFvdw, dFvdxi, JEV] = ...
        WVAMAT(wxi, h, pcs, bcs, joints, Klib);

    Npc = length(pcs);
    Leye = speye(Npts*Nwc);

    Qms = cell(Nh, 1);
    Qvs = cell(Nh, 1);

    dQmdws = cell(Nh, 1);
    dQvdws = cell(Nh, 1);

    dQmdxis = cell(Nh, 1);
    dQvdxis = cell(Nh, 1);

    % Tms = cell(Nh, 1);

    Amatr = cell(Nh, 3);
    Fvr = cell(Nh, 3);
    
    Ljr = cell(length(JEV), 1);
    dLjdwr = cell(length(JEV), 1);
    dLjdxir = cell(length(JEV), 1);
    Ljvr = cell(length(JEV), 1);
    dLjvdwr = cell(length(JEV), 1);
    dLjvdxir = cell(length(JEV), 1);
    Gjr = cell(length(JEV), 1);
    dGjdwr = cell(length(JEV), 1);
    dGjdxir = cell(length(JEV), 1);
    for i=1:length(JEV)
        Ljr{i} = cell(Nh, 1);
        dLjdwr{i} = cell(Nh, 1);
        dLjdxir{i} = cell(Nh, 1);
        Ljvr{i} = cell(Nh, 1);
        dLjvdwr{i} = cell(Nh, 1);
        dLjvdxir{i} = cell(Nh, 1);
        Gjr{i} = cell(Nh, 1);
        dGjdwr{i} = cell(Nh, 1);
        dGjdxir{i} = cell(Nh, 1);
    end

    used_eqis = [];
    for ih=1:Nh
        hstart = (ih-1)*Npts*Nwc;
        k = 0;

        Qm = zeros(Npts*Nwc, Npc*Nwc);
        dQmdw = zeros(Npts*Nwc, Npc*Nwc, Nc);
        dQmdxi = zeros(Npts*Nwc, Npc*Nwc, Npar);
        Qv = zeros(Npts*Nwc, 1);
        dQvdw = zeros(Npts*Nwc, Nc);
        dQvdxi = zeros(Npts*Nwc, Npar);

        % 1. Simplify All Propagations
        for i=1:length(pcs)
            qi = (i-1)*Nwc+1;
            qe = i*Nwc;

            ri = ((pcs(i).irange(1)-1)+1-1)*Nwc+1;
            re = ((pcs(i).irange(1)-1)+1)*Nwc;

            Qm(ri:re, qi:qe) = speye(Nwc);
            for n=1:pcs(i).N-1
                k = k+1;
                si = (k-1)*Nwc+1;
                se = k*Nwc;

                ri = ((pcs(i).irange(1)-1)+n-1)*Nwc+1;
                re = ((pcs(i).irange(1)-1)+n)*Nwc;

                rpi = ((pcs(i).irange(1)-1)+n)*Nwc+1;
                rpe = ((pcs(i).irange(1)-1)+n+1)*Nwc;

                Qm(rpi:rpe, qi:qe) = Amat(hstart+(si:se), hstart+(ri:re))*Qm(ri:re, qi:qe);
                Qv(rpi:rpe) = Amat(hstart+(si:se), hstart+(ri:re))*Qv(ri:re) - ...
                    Fv(hstart+(si:se));

                dQmdw(rpi:rpe, qi:qe, :) = cell2mat( arrayfun(...
                    @(iw) dAmatdw(hstart+(si:se), hstart+(ri:re), iw)*Qm(ri:re, qi:qe) + ...
                    Amat(hstart+(si:se), hstart+(ri:re))*dQmdw(ri:re, qi:qe, iw), ...
                    permute(1:Nc, [1 3 2]), 'UniformOutput', false) );
                dQmdxi(rpi:rpe, qi:qe, :) = cell2mat( arrayfun(...
                    @(ix) dAmatdxi(hstart+(si:se), hstart+(ri:re), ix)*Qm(ri:re, qi:qe) + ...
                    Amat(hstart+(si:se), hstart+(ri:re))*dQmdxi(ri:re, qi:qe, ix), ...
                    permute(1:Npar, [1 3 2]), 'UniformOutput', false) );

                dQvdw(rpi:rpe) = cell2mat( arrayfun( ...
                    @(iw) dAmatdw(hstart+(si:se), hstart+(ri:re), iw)*Qv(ri:re) + ...
                    Amat(hstart+(si:se), hstart+(ri:re))*dQvdw(ri:re, iw) - ...
                    dFvdw(hstart+(si:se), iw), 1:Nc, 'Uniformoutput', false) );
                dQvdxi(rpi:rpe) = cell2mat( arrayfun( ...
                    @(ix) dAmatdxi(hstart+(si:se), hstart+(ri:re), ix)*Qv(ri:re) + ...
                    Amat(hstart+(si:se), hstart+(ri:re))*dQvdxi(ri:re, ix) - ...
                    dFvdw(hstart+(si:se), ix), 1:Npar, 'Uniformoutput', false) );
            end
        end
        ktn = (Npts-length(pcs))*Nwc;
        used_eqis = (1:ktn);
        
        % 2. Simplify All Boundary Conditions & linear joints
        ktn = hstart + (Npts-length(pcs))*Nwc;
        nbcs = sum([bcs.nof]);
        eqis = (1:nbcs);

        % linear joints
        linjis = find(arrayfun(@(j) isempty(j.nl), joints));  % Linear joints

        jeqis = cumsum([1 joints.nof]);  % Joint eq ids
        j0 = jeqis(1:end-1);
        je = jeqis(2:end)-1;

        lis = [];
        for i=linjis
            lis = [lis j0(i):je(i)];
        end
        eqis = [eqis nbcs+lis];
        
        B1 = Amat(ktn+eqis, hstart+(1:Npts*Nwc))*Qm;
        B2 = Amat(ktn+eqis, hstart+(1:Npts*Nwc))*Qv - Fv(ktn+eqis);

        dB1dw = cell2mat( arrayfun( ...
            @(iw) dAmatdw(ktn+eqis, hstart+(1:Npts*Nwc), iw)*...
            Qm(1:Npts*Nwc, 1:Npc*Nwc) + ...
            Amat(ktn+eqis, hstart+(1:Npts*Nwc))*...
            dQmdw(1:Npts*Nwc, 1:Npc*Nwc, iw), ...
            permute(1:Nc, [1 3 2]), 'Uniformoutput', false));
        dB1dxi = cell2mat( arrayfun( ...
            @(ixi) dAmatdxi(ktn+eqis, hstart+(1:Npts*Nwc), ixi)*...
            Qm(1:Npts*Nwc, 1:Npc*Nwc) + ...
            Amat(ktn+eqis, hstart+(1:Npts*Nwc))*...
            dQmdxi(1:Npts*Nwc, 1:Npc*Nwc, ixi), ...
            permute(1:Npar, [1 3 2]), 'Uniformoutput', false) );
        
        dB2dw = cell2mat( arrayfun( ...
            @(iw) dAmatdw(ktn+eqis, hstart+(1:Npts*Nwc),iw)*...
            Qv(1:Npts*Nwc) + Amat(ktn+eqis, hstart+(1:Npts*Nwc))*...
            dQvdw(1:Npts*Nwc,iw) - dFvdw(ktn+eqis,iw), ...
            1:Nc, 'Uniformoutput', false) );
        dB2dxi = cell2mat( arrayfun( ...
            @(ixi) dAmatdxi(ktn+eqis, hstart+(1:Npts*Nwc),ixi)*...
            Qv(1:Npts*Nwc) + Amat(ktn+eqis, hstart+(1:Npts*Nwc))*...
            dQvdxi(1:Npts*Nwc,ixi) - dFvdxi(ktn+eqis,ixi), ...
            1:Nc, 'Uniformoutput', false) );

        % Analytical Null-Space
        neqs = length(eqis);
        rng(1);
        ninds = find(any(B1));
        while rank(B1(:, ninds(1:neqs)))<neqs
            ninds = ninds(randperm(length(ninds)));
        end
        ninds = [ninds setdiff(1:size(B1,2), ninds)];
        [~, si] = sort(ninds);

        tmp1 = B1(:, ninds(1:neqs))\B1(:, ninds(neqs+1:end));
        tmp2 = B1(:, ninds(1:neqs))\B2;

        dtmp1dw = cell2mat( arrayfun( ...
            @(iw) -B1(:, ninds(1:neqs))\dB1dw(:, ninds(1:neqs), iw)*tmp1 + ...
            B1(:, ninds(1:neqs))\dB1dw(:, ninds(neqs+1:end), iw), ...
            permute(1:Nc, [1 3 2]), 'Uniformoutput', false) );
        dtmp1dxi = cell2mat( arrayfun( ...
            @(ixi) -B1(:, ninds(1:neqs))\dB1dxi(:, ninds(1:neqs), ixi)*tmp1 + ...
            B1(:, ninds(1:neqs))\dB1dxi(:, ninds(neqs+1:end), ixi), ...
            permute(1:Npar, [1 3 2]), 'Uniformoutput', false) );

        dtmp2dw = cell2mat( arrayfun( ...
            @(iw) -B1(:, ninds(1:neqs), iw)\dB1dw(:, ninds(1:neqs))*tmp2 + ...
            B1(:, ninds(1:neqs))\dB2dw, 1:Nc, 'Uniformoutput', false) );
        dtmp2dxi = cell2mat( arrayfun( ...
            @(ixi) -B1(:, ninds(1:neqs), ixi)\dB1dxi(:, ninds(1:neqs))*tmp2 + ...
            B1(:, ninds(1:neqs))\dB2dxi, 1:Nc, 'Uniformoutput', false) );
        
        Bm = [-tmp1; speye(size(B1,2)-neqs)];
        Bv = [-tmp2; zeros(size(B1,2)-neqs,1)];
        Bm = Bm(si, :);
        Bv = Bv(si, :);

        dBmdw = cat(1, -dtmp1dw, zeros(size(B1,2)-neqs));
        dBvdw = cat(1, -dtmp2dw, zeros(size(B1,2)-neqs,1));
        dBmdxi = cat(1, -dtmp1dxi, zeros(size(B1,2)-neqs));
        dBvdxi = cat(1, -dtmp2dxi, zeros(size(B1,2)-neqs,1));

        dBmdw = dBmdw(si, :, :);
        dBvdw = dBvdw(si, :);

        dQvdw = cell2mat( arrayfun( ...
            @(iw) dQmdw(:,:,iw)*Bv+Qm*dBvdw(:,iw), permute(1:Nc, [1 3 2]), ...
            'Uniformoutput', false) ) + dQvdw;
        dQvdxi = cell2mat( arrayfun( ...
            @(ixi) dQmdxi(:,:,ixi)*Bv+Qm*dBvdxi(:,ixi), permute(1:Nc, [1 3 2]), ...
            'Uniformoutput', false) ) + dQvdxi;        
        
        dQmdw = cell2mat( arrayfun( ...
            @(iw) dQmdw(:,:,iw)*Bm+Qm*dBmdw(:,:,iw), permute(1:Nc, [1 3 2]), ...
            'Uniformoutput', false) );
        dQmdxi = cell2mat( arrayfun( ...
            @(ixi) dQmdxi(:,:,ixi)*Bm+Qm*dBmdxi(:,:,ixi), permute(1:Nc, [1 3 2]), ...
            'Uniformoutput', false) );
        
        Qv = Qm*Bv+Qv;
        Qm = Qm*Bm;

        used_eqis = [used_eqis (Npts-length(pcs))*Nwc+eqis];
        unused_eqis = setdiff(1:Npts*Nwc, used_eqis);

        %% Finalize
        Qms{ih} = Qm;
        Qvs{ih} = Qv;

        dQmdws{ih} = dQmdw;
        dQmdxis{ih} = dQmdxi;

        dQvdws{ih} = dQvdw;
        dQvdxis{ih} = dQvdxi;

        %% Transform
        hisf = hstart+(1:Npts*Nwc);
        his = hstart+unused_eqis;
        Amatr{ih, 1} = Amat(his,hisf)*Qms{ih};
        Amatr{ih, 2} = dAmatdw(his,hisf)*Qms{ih}+Amat(his,hisf)*dQmdws{ih};
        Amatr{ih, 3} = dAmatdxi(his,hisf)*Qms{ih}+Amat(his,hisf)*dQmdxis{ih};
        
        Fvr{ih, 1} = Fv(his) - Amat(his,hisf)*Qvs{ih};
        Fvr{ih, 2} = dFvdw(his)-dAmatdw(his,hisf)*Qvs{ih}-Amat(his,hisf)*dQvdws{ih};
        Fvr{ih, 3} = dFvdxi(his)-dAmatdxi(his,hisf)*Qvs{ih}-Amat(his,hisf)*dQvdxis{ih};

        for i=1:length(JEV)
            nd = size(JEV.Lj, 1)/Nh;
            hisn = (ih-1)*nd+(1:nd);
            
            Ljr{i}{ih} = JEV.Lj(hisn,hisf)*Qms{ih};
            dLjdwr{i}{ih} = JEV.dLjdw(hisn,hisf)*Qms{ih} + JEV.Lj(hisn,hisf)*dQmdws{ih};
            dLjdxir{i}{ih} = JEV.dLjdxi(hisn,hisf)*Qms{ih} + JEV.Lj(hisn,hisf)*dQmdxis{ih};

            Ljvr{i}{ih} = JEV.Lj(hisn,hisf)*Qvs{ih};
            dLjvdwr{i}{ih} = JEV.dLjdw(hisn,hisf)*Qvs{ih} + JEV.Lj(hisn,hisf)*dQvdws{ih};
            dLjvdxir{i}{ih} = JEV.dLjdxi(hisn,hisf)*Qvs{ih} + JEV.Lj(hisn,hisf)*dQvdxis{ih};

            nf = size(JEV.Gj, 2)/Nh;
            hisn = (ih-1)*nf+(1:nf);
            
            Gjr{i}{ih} = JEV.Gj(his, hisn);
            dGjdwr{i}{ih} = JEV.dGjdw(his, hisn);
            dGjdxir{i}{ih} = JEV.dGjdxi(his, hisn);
        end

        if any(h(ih,:)~=0) && ~isempty(varargin)
            if strcmp(varargin{1}, 'r')
                tmp = cat(2, Qms{ih}, 1j*Qms{ih});
                Qms{ih} = cat(1, real(tmp), imag(tmp));

                tmp = cat(2, dQmdws{ih}, 1j*dQmdws{ih});
                dQmdws{ih} = cat(1, real(tmp), imag(tmp));

                tmp = cat(2, dQmdxis{ih}, 1j*dQmdxis{ih});
                dQmdxis{ih} = cat(1, real(tmp), imag(tmp));

                Qvs{ih} = cat(1, real(Qvs{ih}), imag(Qvs{ih}));
                dQvdws{ih} = cat(1, real(dQvdws{ih}), imag(dQvdws{ih}));
                dQvdxis{ih} = cat(1, real(dQvdxis{ih}), imag(dQvdxis{ih}));
            
                for i=1:3
                    tmp = cat(2, Amatr{ih, i}, 1j*Amatr{ih, i});
                    Amatr{ih, i} = cat(1, real(tmp), imag(tmp));
                    
                    Fvr{ih, i} = cat(1, real(Fvr{ih, i}), imag(Fvr{ih, i}));
                end
                for i=1:length(JEV)
                    tmp = cat(2, Ljr{i}{ih}, 1j*Ljr{i}{ih});
                    Ljr{i}{ih} = cat(1, real(tmp), imag(tmp));
                    
                    tmp = cat(2, dLjdwr{i}{ih}, 1j*dLjdwr{i}{ih});
                    dLjdwr{i}{ih} = cat(1, real(tmp), imag(tmp));
                    
                    tmp = cat(2, dLjdxir{i}{ih}, 1j*dLjdxir{i}{ih});
                    dLjdxir{i}{ih} = cat(1, real(tmp), imag(tmp));
                    
                    Ljvr{i}{ih} = cat(1, real(Ljvr{i}{ih}), imag(Ljvr{i}{ih}));
                    dLjvdwr{i}{ih} = cat(1, real(dLjvdwr{i}{ih}), imag(dLjvdwr{i}{ih}));
                    dLjvdxir{i}{ih} = cat(1, real(dLjvdxir{i}{ih}), imag(dLjvdxir{i}{ih}));
                    
                    tmp = cat(2, Gjr{i}{ih}, 1j*Gjr{i}{ih});
                    Gjr{i}{ih} = cat(1, real(tmp), imag(tmp));
                    
                    tmp = cat(2, dGjdwr{i}{ih}, 1j*dGjdwr{i}{ih});
                    dGjdwr{i}{ih} = cat(1, real(tmp), imag(tmp));
                    
                    tmp = cat(2, dGjdxir{i}{ih}, 1j*dGjdxir{i}{ih});
                    dGjdxir{i}{ih} = cat(1, real(tmp), imag(tmp));
                end
            end
        end            
    end
    
    dAmatdxir = blkdiag(Amatr{:, 3});
    dAmatdwr = blkdiag(Amatr{:, 2});
    Amatr = blkdiag(Amatr{:, 1});

    dFvdxir = cell2mat(Fvr(:, 3));
    dFvdwr = cell2mat(Fvr(:, 2));
    Fvr = cell2mat(Fvr(:, 1));

    for i=1:length(JEV)
        Ljr{i} = blkdiag(Ljr{i}{:});
        dLjdwr{i} = blkdiag(dLjdwr{i}{:});
        dLjdxir{i} = blkdiag(dLjdxir{i}{:});

        Ljvr{i} = cell2mat(Ljvr{i});
        dLjvdwr{i} = cell2mat(dLjvdwr{i});
        dLjvdxir{i} = cell2mat(dLjvdxir{i});

        Gjr{i} = blkdiag(Gjr{i}{:});
        dGjdwr{i} = blkdiag(dGjdwr{i}{:});
        dGjdxir{i} = blkdiag(dGjdxir{i}{:});
    end    

    JEVr = struct('Lj', Ljr, 'Ljv', Ljvr, 'dLjdw', dLjdwr, 'dLjdxi', dLjdxir, ...
                  'dLjvdw', dLjvdwr, 'dLjvdxi', dLjvdxir, 'Gj', Gjr, 'dGjdw', dGjdwr, ...
                  'dGjdxi', dGjdxir);

    RECOV = struct('Qm', blkdiag(Qms{:}), 'dQmdw', blkdiag(dQmdws{:}), ...
                   'dQmdxi', blkdiag(dQmdxis{:}), 'Qv', cell2mat(Qvs), ...
                   'dQvdw', cell2mat(dQvdws), 'dQvdxi', cell2mat(dQvdxis));
end
