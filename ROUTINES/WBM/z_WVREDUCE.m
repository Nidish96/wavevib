function [Qreds, dQredsdw, dQredsdxi, Treds] = WVREDUCE(wxi, h, pcs, bcs, joints, Klib)
%WVREDUCE Returns reduction relationships for a WBM.
%
%	USAGE:
%		[Qred, Tred] = WVREDUCE(wxi, h, pcs, bcs, joints, Klib);
%	INPUTS:
%		wxi, h, pcs, bcs, joints, Klib
%	OUTPUTS:
%		Qred, Tred

    Nwc = size(pcs(1).wcomps,1);
    assert(all(arrayfun(@(p) size(p.wcomps,1), pcs(2:end))==Nwc));
    Npts = pcs(end).irange(end);
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
    Nh = size(h, 1);
    Nc = size(h, 2);
    Npar = length(wxi)-Nc;
    assert(Npar==1);
    
    % Nonlinearities
    nlis = find(arrayfun(@(j) ~isempty(j.nl), joints));
    nnl = length(nlis);
    
    if nnl==0
        % Model is linear
        Qred = nan;
        Tred = nan;
        return
    end

    bcis = [bcs.i];
    Leye = speye(Npts);

    % Model has nonlinear joints.
    % Full model can be expressed in terms of joints
    redset = struct('is', [], 'ps', []);
    bcjs = cell(1,nnl);
    for i=1:nnl
        k = nlis(i);

        redset.is = [redset.is joints(k).is(1:joints(k).type)];
        redset.ps = [redset.ps joints(k).ps(1:joints(k).type)];

        bcjs{i} = find(any(bcis(:)==joints(k).is(:)',2));  % bcs involving joint pts
    end
    assert(length(unique(redset.ps))==length(redset.ps));
    % one piece can't have multiple joints
    Nptred = length(redset.is);

    Qreds = cell(Nh, 1);
    dQredsdxi = cell(Nh, 1);
    dQredsdw = cell(Nh, 1);

    %% Zeroth Harmonic
    i1 = 1;
    if all(h(1,:)==0) % Zeroth Harmonic 
        rmats = {};
        drmatsdxi = {};
        for i=1:nnl
            k = nlis(i);
            
            if isa(joints(k).nlfcofs0, 'function_handle')
                nlfcfs = joints(k).nlfcofs0(wxi(Nc+1:end));
            else
                nlfcfs = joints(k).nlfcofs0;
            end
            alzi = find(all(nlfcfs==0,2));
            nzi = length(alzi);

            if isa(joints(k).cofs0, 'function_handle')
                jcofs = joints(k).cofs0(wxi(Nc+1:end));
                djcofsdxi = joints(k).dcofsdxi0(wxi(Nc+1:end));
            else
                jcofs = joints(k).cofs0;
                djcofsdxi = zeros(size(jcofs));
            end

            bmats = [];
            dbmatsdxi = [];
            for bi=1:length(bcjs{i})
                if isa(bcs(bcjs{i}(bi)).cofs0, 'function_handle')
                    bmats = [bmats;
                             kron(Leye(bcs(bcjs{i}(bi)).i, joints(k).is), ...
                                  bcs(bcjs{i}(bi)).cofs0(wxi(1:Nc), wxi(Nc+1:end)))];
                    dbmatsdxi = [dbmatsdxi;
                                 kron(Leye(bcs(bcjs{i}(bi)).i, joints(k).is), ...
                                      bcs(bcjs{i}(bi)).dcofsdxi0(wxi(1:Nc), wxi(Nc+1:end)))];
                else
                    tmp = kron(Leye(bcs(bcjs{i}(bi)).i, joints(k).is), bcs(bcjs{i}(bi)).cofs0);
                    bmats = [bmats; tmp];
                    dbmatsdxi = [dbmatsdxi; zeros(size(tmp))];
                end
            end
            jmat = full([jcofs(alzi,:); bmats]);
            djmatdxi = [djcofsdxi(alzi,:); dbmatsdxi];

            ninds = find(any(jmat));
            rng(1)
            while rank(jmat(:, ninds(1:Nwc*nzi)))<Nwc*nzi
                ninds = ninds(randperm(length(ninds)));
            end
            ninds = [ninds setdiff(1:size(jmat,2), ninds)];
            [~, si] = sort(ninds);

            tmp = jmat(:, ninds(1:nzi*Nwc))\jmat(:, ninds(nzi*Nwc+1:end));

            dtmpdxi = -(jmat(:, ninds(1:nzi*Nwc))\djmatdxi(:, ninds(1:nzi*Nwc)))*tmp + ...
                      jmat(:, ninds(1:nzi*Nwc))\djmatdxi(:, ninds(nzi*Nwc+1:end));
            nlJ = [-tmp; kron(speye(Nwc), speye(joints(k).type-nzi))];
            dnlJdxi = [-dtmpdxi; sparse(zeros((joints(k).type-nzi)*Nwc))];

            rmats = cat(1, rmats, full(nlJ(si, :)));
            drmatsdxi = cat(1, drmatsdxi, full(dnlJdxi(si, :)));
        end

        Qred = zeros(Npts*Nwc, Nptred*Nwc);
        dQreddxi = zeros(Npts*Nwc, Nptred*Nwc);

        for i=1:Nptred
            m = redset.is(i);
            qi = (i-1)*Nwc+1;
            qe = i*Nwc;
            
            pci = redset.ps(i);
            ml = m-pcs(pci).irange(1)+1;

            for n=pcs(pci).irange(1):pcs(pci).irange(end)
                si = (n-1)*Nwc+1;
                se = n*Nwc;

                nl = n-pcs(pci).irange(1)+1;

                ns = (0:Nwc-1);
                ks = (0:Nwc-1)';

                dx = sign(ml-nl+eps)*abs(diff(pcs(pci).U([nl ml]))*pcs(pci).S);
                AA = triu((-1).^(ns-ks).*factorial(ns)./(factorial(max(ns-ks,0)).*factorial(ks)).*(dx.^(ns-ks)));

                Qred(si:se, qi:qe) = AA;
            end
        end
        Qreds{i1} = Qred*blkdiag(rmats{:});
        dQredsdxi{i1} = Qred*blkdiag(drmatsdxi{:});
        dQredsdw{i1} = zeros(size(Qreds{i1}));
        
        i1 = 2;
    end

    %% Non-Zero Harmonics
    
    for ih=i1:Nh
        rmats = {};
        drmatsdw = {};
        drmatsdxi = {};
        for i=1:nnl
            k = nlis(i);
        
            if isa(joints(k).nlfcofs, 'function_handle')
                nlfcfs = joints(k).nlfcofs(h(ih,:)*wxi(1:Nc), wxi(Nc+1:end));
            else
                nlfcfs = joints(k).nlfcofs;
            end
            alzi = find(all(nlfcfs==0,2));
            nzi = length(alzi);

            if isa(joints(k).cofs, 'function_handle')
                jcofs = joints(k).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+1:end));
                djcofsdw = joints(k).dcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+1:end)).*...
                    permute(h(ih,:), [1 3 2]);
                djcofsdxi = joints(k).dcofsdxi(h(ih,:)*wxi(1:Nc), wxi(Nc+1:end));
            else
                jcofs = joints(k).cofs;
                djcofsdw = zeros(size(jcofs,1), size(jcofs,2), Nc);
                djcofsdxi = zeros(size(jcofs,1), size(jcofs,2));
            end

            bmats = [];
            dbmatsdw = [];
            dbmatsdxi = [];
            for bi=1:length(bcjs{i})
                if isa(bcs(bcjs{i}(bi)).cofs, 'function_handle')
                    bmats = [bmats;
                             kron(Leye(bcs(bcjs{i}(bi)).i, joints(k).is), ...
                                  bcs(bcjs{i}(bi)).cofs(h(ih,:)*wxi(1:Nc), wxi(Nc+1:end)))];
                    dbmatsdxi = [dbmatsdxi;
                                 kron(Leye(bcs(bcjs{i}(bi)).i, joints(k).is), ...
                                      bcs(bcjs{i}(bi)).dcofsdxi0(wxi(1:Nc), wxi(Nc+1:end)))];
                    dbmatsdw = cat(3, dbmatsdw, ...
                                   kron(Leye(bcs(bcjs{i}(bi)).i, joints(k).is), ...
                                        bcs(bcjs{i}(bi)).dcofsdw(h(ih,:)*wxi(1:Nc), wxi(Nc+1:end))).*permute(h(ih,:), [1 3 2]));
                else
                    tmp = kron(Leye(bcs(bcjs{i}(bi)).i, joints(k).is), bcs(bcjs{i}(bi)).cofs);
                    bmats = [bmats; tmp];
                    dbmatsdxi = [dbmatsdxi; zeros(size(tmp))];
                    dbmatsdw = cat(1, dbmatsdw, zeros(size(tmp,1), size(tmp,2), Nc));
                end
            end
            jmat = [jcofs(alzi,:); bmats];
            djmatdw = cat(1, djcofsdw(alzi,:,:), dbmatsdw);
            djmatdxi = [djcofsdxi(alzi,:); dbmatsdxi];
            
            ninds = ((1:joints(k).type)'-1)*Nwc+(1:Nwc);

            tmp = jmat(:, ninds(1:nzi*Nwc))\jmat(:, ninds(nzi*Nwc+1:end));
            
            dtmpdw = cell2mat(arrayfun(@(ia) jmat(:, ninds(1:nzi*Nwc))\...
                                       (djmatdw(:, ninds(1:nzi*Nwc), ia)*tmp + ...
                                        djmatdw(:, ninds(nzi*Nwc+1:end))), ...
                                       permute(1:Nc, [1 3 2]), 'UniformOutput', false));
            dtmpdxi = -jmat(:, ninds(1:nzi*Nwc))\(djmatdxi(:, ninds(1:nzi*Nwc))*tmp + ...
                                                  djmatdxi(:, ninds(nzi*Nwc+1:end)));
            nlJ = zeros(size(tmp)+[Nwc*(joints(k).type-nzi) 0]);
            dnlJdw = zeros(size(tmp,1)+Nwc*(joints(k).type-nzi), size(tmp,2), Nc);
            dnlJdxi = zeros(size(tmp)+[Nwc*(joints(k).type-nzi) 0]);
            nlJ([ninds(nzi*Nwc+1:end) ninds(1:nzi*Nwc)], :) = ...
                [kron(speye(Nwc), speye(joints(k).type-nzi)); -tmp];
            dnlJdw([ninds(nzi*Nwc+1:end) ninds(1:nzi*Nwc)], :, :) = ...
                cat(1, zeros(Nwc*(joints(k).type-nzi), Nwc*(joints(k).type-nzi), Nc), ...
                    -dtmpdw);

            rmats = cat(1, rmats, nlJ);
            drmatsdxi = cat(1, drmatsdxi, dnlJdxi);
            dnlJdw = mat2cell(dnlJdw, size(dnlJdw,1), size(dnlJdw,2), ones(1,Nc));
            drmatsdw = cat(1, drmatsdw, dnlJdw);
        end
        
        Qred = zeros(Npts*Nwc, Nptred*Nwc);
        dQreddw = zeros(Npts*Nwc, Nptred*Nwc, Nc);
        dQreddxi = zeros(Npts*Nwc, Nptred*Nwc);
        for i=1:Nptred
            m = redset.is(i);
            qi = (i-1)*Nwc+1;
            qe = i*Nwc;
            
            pci = redset.ps(i);
            ml = m-pcs(pci).irange(1)+1;
            
            Ks = arrayfun(@(k) k.K(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))), Klib);  %(nk,1)
            dKdws = cell2mat(arrayfun(@(k) k.dKdw(h(ih,:)*wxi(1:Nc), ...
                                                  wxi(Nc+(1:Npar)))*h(ih, :), ...
                                      Klib, 'UniformOutput', false));  %(nk,Nc)
            dKdxis = arrayfun(@(k) k.dKdxi(h(ih,:)*wxi(1:Nc), ...
                                           wxi(Nc+(1:Npar))), Klib);  %(nk,Npar)

            K = pcs(i).wcomps(:,1).*Ks(pcs(i).wcomps(:,2));  %(Nwc,1)
            dKdw = pcs(i).wcomps(:,1).*dKdws(pcs(i).wcomps(:,2),:);  %(Nwc,Nc)
            dKdxi = pcs(i).wcomps(:,1).*dKdxis(pcs(i).wcomps(:,2),:);  %(Nwc,Npar)

            for n=pcs(pci).irange(1):pcs(pci).irange(end)
                si = (n-1)*Nwc+1;
                se = n*Nwc;

                nl = n-pcs(pci).irange(1)+1;
                
                dx = sign(ml-nl+eps)*abs(diff(pcs(pci).U([nl ml]))*pcs(pci).S);

                Qred(si:se, qi:qe) = diag(exp(-K*dx));
                dQreddw(si:se, qi:qe, :) = reshape(cell2mat(arrayfun(@(a) -diag(exp(-K*dx)*dx.*dKdw(:,a)), 1:Nc, 'UniformOutput', false)), Nwc,Nwc,Nc);
                dQreddxi(si:se, qi:qe, :) = reshape(cell2mat(arrayfun(@(a) -diag(exp(-K*dx)*dx.*dKdxi(:,a)), 1:Npar, 'UniformOutput', false)), Nwc,Nwc,Npar);
            end
        end
        rmats = blkdiag(rmats{:});
        Qredr = Qred*rmats;
        
        dQreddw = cell2mat(arrayfun(@(iw) dQreddw(:,:,iw)*rmats + ...
                                    Qred*blkdiag(drmatsdw{:, iw}), ...
                                    permute(1:Nc, [1 3 2]), 'Uniformoutput', false));
        dQreddxi = dQreddxi*rmats + Qred*blkdiag(drmatsdxi{:});

        tmp = [Qredr 1j*Qredr];
        Qreds{ih} = [real(tmp); imag(tmp)];

        tmp = cat(2, dQreddw, 1j*dQreddw);
        dQredsdw{ih} = cat(1, real(tmp), imag(tmp));

        tmp = [dQreddxi 1j*dQreddxi];
        dQredsdxi{ih} = [real(tmp); imag(tmp)];
    end

    Qreds = blkdiag(Qreds{:});
    dQredsdw = blkdiag(dQredsdw{:});
    dQredsdxi = blkdiag(dQredsdxi{:});

    % %% Accounting for Excitation
    % if all(h(1,:)==0)
    %     Qv = mat2cell(zeros(Npts*Nwc*Nhc,1), [Npts*Nwc, repmat(2*Npts*Nwc,1,Nh-1)], 1);
    %     dQvdxi = mat2cell(zeros(Npts*Nwc*Nhc,Npar), [Npts*Nwc, repmat(2*Npts*Nwc,1,Nh-1)], Npar);
    %     dQvdw = mat2cell(zeros(Npts*Nwc*Nhc,Nc), [Npts*Nwc, repmat(2*Npts*Nwc,1,Nh-1)], Nc);
    % else
    %     Qv = mat2cell(zeros(Npts*Nwc*Nhc,1), [repmat(2*Npts*Nwc,1,Nh)], 1);
    %     dQvdxi = mat2cell(zeros(Npts*Nwc*Nhc,Npar), [repmat(2*Npts*Nwc,1,Nh)], Npar);
    %     dQvdw = mat2cell(zeros(Npts*Nwc*Nhc,Nc), [repmat(2*Npts*Nwc,1,Nh)], Nc);
    % end

    % k = 0;
    % for i=1:length(pcs)  % Excitation points
    %     for n=1:pcs(i).N-1
    %         k = k+1;
    %         si = (k-1)*Nwc+1;
    %         se = k*Nwc;
            
    %         dx = abs(diff(pcs(i).U(n:n+1))*pcs(i).S);
    %         ei = find(pcs(i).exci==n);
    %         if dx==0 && ~isempty(ei)
    %             assert(length(ei)==1);
                
    %             if all(pcs(i).excnh(ei,:)==0)  % zeroth harmonic excitation
    %                 if numel(pcs(i).exccofs0{ei})==1
    %                     Qv{1}(si:se) = pcs(i).exccofs0{ei};
    %                 else
    %                     Qv{1}(si:se) = pcs(i).exccofs0{ei}{1}(wxi(Nc+(1:Npar)));
    %                     dQvdxi{1}(si:se,:) = pcs(i).exccofs0{ei}{1}(wxi(Nc+(1:Npar)));
    %                 end
    %             else
    %                 ih = find(all(pcs(i).excnh(ei,:)==h, 2));
    %                 if isempty(ih)
    %                     continue;
    %                 end
    %                 if numel(pcs(i).exccofs{ei})==1
    %                     Qv{ih}(si:se) = pcs(i).exccofs{ei};
    %                 else
    %                     tmp = pcs(i).exccofs{ei}{1}(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));
    %                     Qv{ih}([si:se Npts*Nwc+(si:se)]) = [real(tmp);imag(tmp)];
    %                     tmp = pcs(i).exccofs{ei}{2}(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar))).*h(ih,:);

    %                     dQvdw{ih}([si:se Npts*Nwc+(si:se)], :) = [real(tmp); imag(tmp)];
                        
    %                     tmp = pcs(i).exccofs{ei}{3}(h(ih,:)*wxi(1:Nc), wxi(Nc+(1:Npar)));
    %                     dQvdxi{ih}([si:se Npts*Nwc+(si:se)], :) = [real(tmp); imag(tmp)];
    %                 end
    %             end
    %         end
    %     end
    % end
    % Qv = cell2mat(Qv);
    % dQvdw = cell2mat(dQvdw);
    % dQvdxi = cell2mat(dQvdxi);
    
    %% RHS
    bcjs = cell2mat(bcjs);
    bcjs_ = setdiff(1:length(bcs), bcjs);

    ktn = (Npts-length(pcs))*Nwc;

    bnof = cumsum([1 bcs.nof]);
    b0 = bnof(1:end-1);
    be = bnof(2:end)-1;

    Leye = speye(Npts*Nwc);
    Tred = [];
    for i=1:length(bcjs_)
        Tred = [Tred; Leye(ktn+(b0(bcjs_(i)):be(bcjs_(i))),:)];
    end
    ktn = ktn+sum([bcs.nof]);
    for i=1:nnl
        k = nlis(i);
        
        Tred = [Tred; Leye(ktn+(1:joints(k).type-1), :)];
    end

    Treds = kron(speye(Nhc), Tred);
end
