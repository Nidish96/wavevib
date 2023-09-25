function [pcs, bcs, joints, excs, Klib] = WBPREPROC(pcs, bcs, joints, excs, Klib)

    origcoords = cell2mat({pcs.coords}');
    if size(origcoords,2)<3
        origcoords = [origcoords zeros(size(origcoords,1),3)];
        origcoords = origcoords(:,1:3);
    end
    origids = cell2mat(arrayfun(@(i) ones(size(pcs(i).coords,1),1)*i, 1:length(pcs), 'UniformOutput', false)');
    origcoords = [origcoords origids];  % last column is piece ID (1-indexed)
    syms w xi
    % 0. Preprocess Klib (for derivatives)
    for i=1:length(Klib)
        if ~isfield(Klib(i), 'dKdw') || isempty(Klib(i).dKdw)
            Klib(i).dKdw = matlabFunction(diff(Klib(i).K(w,xi), w), 'Vars', [w xi]);
        end
        if ~isfield(Klib(i), 'dKdxi') || isempty(Klib(i).dKdxi)
            Klib(i).dKdxi = matlabFunction(diff(Klib(i).K(w,xi), xi), 'Vars', [w xi]);
        end
    end
    % 1. Preprocess Pieces
    Nbef = 0;
    elist = {'exci', 'excnh'};
    for i=1:length(pcs)
        if size(pcs(i).coords,3)<3
            pcs(i).coords = [pcs(i).coords zeros(size(pcs(i).coords,1),3)];
            pcs(i).coords = pcs(i).coords(:, 1:3);
        end
        coords0 = pcs(i).coords(1,:);
        [U, S, V] = svd(pcs(i).coords-coords0,0);
        if nnz(S)>1  % rank(pcs(i).coords)>1
            pcs(i).err = true;
            pcs(i).msg = sprintf('Points on Piece %d are not on a line. Check the coordinates.',i);
            warning('Points on Piece %d are not on a line. Check the coordinates.',i);
            continue;
        end
        pcs(i).err = false;
        U = U(:,1); S = S(1); V = V(:,1);
        if any(U~=sort(U)) && any(U~=sort(U,'descend'))
            pcs(i).err = true;
            pcs(i).msg = sprintf('Points on Piece %d are not in sequence. Check the coordinates.', i);
            warning('Points on Piece %d are not in sequence. Check the coordinates.', i);
        end
        if length(U)~=((length(unique(U))-2)*2+2)
            U = U([1; kron((2:length(U)-1)',[1;1]); end]);
            pcs(i).coords = U*S*V'+coords0;
        end
        for ei=elist
            if ~isfield(pcs(i), ei{:})
                pcs(i).(ei{:}) = [];
            end
        end
        if ~isfield(pcs(i), 'exccofs') || ~iscell(pcs(i).exccofs)
            pcs(i).exccofs = {};
        end
        if ~isfield(pcs(i), 'exccofs0') || ~iscell(pcs(i).exccofs0)
            pcs(i).exccofs0 = {};
        end

        pcs(i).N = size(pcs(i).coords,1);  % Number of Elements
        pcs(i).L = vecnorm(diff(pcs(i).coords([1 end],:)));
        pcs(i).nhat = diff(pcs(i).coords([1 end],:))/pcs(i).L;
        pcs(i).S = S*vecnorm(U)*vecnorm(V);
        pcs(i).V = V/vecnorm(V)*sign(U(end)-U(1));
        pcs(i).U = U/vecnorm(U)*sign(U(end)-U(1));

        pcs(i).theta = atan2(vecnorm(pcs(i).nhat(1:2)), pcs(i).nhat(3));
        pcs(i).phi = atan2(pcs(i).nhat(2), pcs(i).nhat(1));
        pcs(i).irange = Nbef + [1 pcs(i).N];
        Nbef = Nbef + pcs(i).N;
    end
    procids = cell2mat(arrayfun(@(i) ones(size(pcs(i).coords,1),1)*i, 1:length(pcs), 'UniformOutput', false)');
    proccoords = [cell2mat({pcs.coords}') procids];

    % 2. Preprocess Boundary Conditions
    Nwc = size(pcs(1).wcomps,1);
    for n=1:length(bcs)
        if ~isfield(bcs(n), 'dcofsdw') || isempty(bcs(n).dcofsdw)
            bcs(n).dcofsdw = matlabFunction(diff(bcs(n).cofs(w,xi), w), 'Vars', [w xi]);
        end
        if ~isfield(bcs(n), 'dcofsdxi') || isempty(bcs(n).dcofsdxi)
            bcs(n).dcofsdxi = matlabFunction(diff(bcs(n).cofs(w,xi), xi), 'Vars', [w xi]);
        end
        bcs(n).nof = size(bcs(n).cofs(w,xi),1);
        if ~isfield(bcs(n), 'rih')
            bcs(n).rih = [];
        end
        if ~isfield(bcs(n), 'rhs') || isempty(bcs(n).rhs)
            bcs(n).rhs = @(w,xi) 0;
            bcs(n).drhsdw = @(w,xi) 0;
            bcs(n).drhsdxi = @(w,xi) 0;
        else
            if ~isfield(bcs(n), 'drhsdw') || isempty(bcs(n).drhsdw)
                bcs(n).drhsdw = matlabFunction(diff(bcs(n).rhs(w,xi), w), 'Vars', [w xi]);
            end
            if ~isfield(bcs(n), 'drhsdxi') || isempty(bcs(n).drhsdxi)
                bcs(n).drhsdxi = matlabFunction(diff(bcs(n).rhs(w,xi), xi), 'Vars', [w xi]);
            end
        end
        
        % Zeroth coefficients
        if ~isfield(bcs(n), 'cofs0') || isempty(bcs(n).cofs0)
            bcs(n).cofs0 = matlabFunction(bcs(n).cofs(0, xi), 'Vars', xi);  % Default behavior
        end
        if ~isfield(bcs(n), 'dcofsdxi0') || isempty(bcs(n).dcofsdxi0)
            bcs(n).dcofsdxi0 = matlabFunction(diff(bcs(n).cofs0(xi), xi), 'Vars', xi);
        end
        assert(size(bcs(n).cofs(w,xi),1) == size(bcs(n).cofs0(xi),1))
        if ~isfield(bcs(n), 'rhs0') || isempty(bcs(n).rhs0)
            bcs(n).rhs0 = @(xi) 0;
            bcs(n).drhsdxi0 = @(xi) 0;
        else
            if ~isfield(bcs(n), 'drhsdxi0') || isempty(bcs(n).drhsdxi0)
                bcs(n).drhsdxi0 = matlabFunction(diff(bcs(n).rhs0(xi), xi), 'Vars', xi);
            end
        end

        % Account for new nodes included in model
        mi = find(all(proccoords==origcoords(bcs(n).i,:), 2));
        bcs(n).i = mi(1);
        bcs(n).pi = find(arrayfun(@(p) prod(bcs(n).i-p.irange)<=0, pcs));
    end

    % 3. Preprocess Joints
    for n=1:length(joints)
        if ~isfield(joints(n), 'cofs') || isempty(joints(n).cofs)
            switch joints(n).type
              case {1, 2}
                joints(n).cofs = @(w,xi) [eye(Nwc) -eye(Nwc)];
              otherwise
                joints(n).cofs = @(w,xi) [eye(Nwc)/joints(n).type -repmat(eye(Nwc), 1,joints(n).type-1)];
            end
        end
        if ~isfield(joints(n), 'dcofsdw') || isempty(joints(n).dcofsdw)
            joints(n).dcofsdw = matlabFunction(diff(joints(n).cofs(w,xi), w), 'Vars', [w xi]);
        end
        if ~isfield(joints(n), 'dcofsdxi') || isempty(joints(n).dcofsdxi)
            joints(n).dcofsdxi = matlabFunction(diff(joints(n).cofs(w,xi), xi), 'Vars', [w xi]);
        end
        joints(n).nof = size(joints(n).cofs(w,xi),1);

        if ~isfield(joints(n), 'rih')
            joints(n).rih = [];
        end
        if ~isfield(joints(n), 'rhs') || isempty(joints(n).rhs)
            joints(n).rhs = @(w,xi) zeros(Nwc,1);
            joints(n).drhsdw = @(w,xi) zeros(Nwc,1);
            joints(n).drhsdxi = @(w,xi) zeros(Nwc,1);
        else
            if ~isfield(joints(n), 'drhsdw') || isempty(joints(n).drhsdw)
                joints(n).drhsdw = matlabFunction(diff(joints(n).rhs(w,xi), w), 'Vars', [w xi]);
            end
            if ~isfield(joints(n), 'drhsdxi') || isempty(joints(n).drhsdxi)
                joints(n).drhsdxi = matlabFunction(diff(joints(n).rhs(w,xi), xi), 'Vars', [w xi]);
            end
        end
        %% Zeroth Harmonics
        if ~isfield(joints(n), 'cofs0') || isempty(joints(n).cofs0)
            joints(n).cofs0 = matlabFunction(joints(n).cofs(0,xi), 'Vars', xi);
        end
        if ~isfield(joints(n), 'dcofsdxi0') || isempty(joints(n).dcofsdxi0)
            joints(n).dcofsdxi0 = matlabFunction(diff(joints(n).cofs0(xi), xi), 'Vars', xi);
        end
        assert(size(joints(n).cofs(w,xi),1) == size(joints(n).cofs0(xi),1));

        if ~isfield(joints(n), 'rhs0') || isempty(joints(n).rhs0)
            joints(n).rhs0 = @(xi) zeros(Nwc,1);
            joints(n).drhsdxi0 = @(xi) zeros(Nwc,1);
        else
            if ~isfield(joints(n), 'drhsdxi0') || isempty(joints(n).drhsdxi0)
                joints(n).drhsdxi0 = matlabFunction(diff(joints(n).rhs0(xi), xi), 'Vars', xi);
            end
        end
        %% Nonlinear Joint
        if ~isfield(joints(n), 'nl') || isempty(joints(n).nl)
            joints(n).nl = [];
        else
            if ~isfield(joints(n), 'dnlfcofsdw') || isempty(joints(n).dnlfcofsdw)
                joints(n).dnlfcofsdw = matlabFunction(diff(joints(n).nlfcofs(w,xi), w), 'Vars', [w xi]);
            end
            if ~isfield(joints(n), 'dnlfcofsdxi') || isempty(joints(n).dnlfcofsdxi)
                joints(n).dnlfcofsdxi = matlabFunction(diff(joints(n).nlfcofs(w,xi), xi), 'Vars', [w xi]);
            end
            
            if ~isfield(joints(n), 'dnldcofsdw') || isempty(joints(n).dnldcofsdw)
                joints(n).dnldcofsdw = matlabFunction(diff(joints(n).nldcofs(w,xi), w), 'Vars', [w xi]);
            end
            if ~isfield(joints(n), 'dnldcofsdxi') || isempty(joints(n).dnldcofsdxi)
                joints(n).dnldcofsdxi = matlabFunction(diff(joints(n).nldcofs(w,xi), xi), 'Vars', [w xi]);
            end

            joints(n).nld = size(joints(n).nldcofs(0,0), 1);

            % Zeroth Harmonic
            if ~isfield(joints(n), 'nlfcofs0') || isempty(joints(n).nlfcofs0)
                joints(n).nlfcofs0 = matlabFunction(joints(n).nlfcofs(0, xi), 'Vars', xi);
            end
            if ~isfield(joints(n), 'dnlfcofsdxi0') || isempty(joints(n).dnlfcofsdxi0)
                joints(n).dnlfcofsdxi0 = matlabFunction(diff(joints(n).nlfcofs0(xi), xi), ...
                                                        'Vars', xi);
            end
            if ~isfield(joints(n), 'nldcofs0') || isempty(joints(n).nldcofs0)
                joints(n).nldcofs0 = matlabFunction(joints(n).nldcofs(0, xi), 'Vars', xi);
            end
            if ~isfield(joints(n), 'dnldcofsdxi0') || isempty(joints(n).dnldcofsdxi0)
                joints(n).dnldcofsdxi0 = matlabFunction(diff(joints(n).nldcofs0(xi), xi), ...
                                                        'Vars', xi);
            end

            assert(size(joints(n).nldcofs(0,0), 1) == size(joints(n).nldcofs0(0), 1));
        end        

        %% Account for new nodes included in model
        switch joints(n).type
          case {1, 2}
            mi = find(all(proccoords==origcoords(joints(n).i,:), 2));
            mj = find(all(proccoords==origcoords(joints(n).j,:), 2));

            joints(n).i = mi(1);
            mj = setdiff(mj, mi(1));
            joints(n).j = mj(1);

            joints(n).pi = find(arrayfun(@(p) prod(joints(n).i-p.irange)<=0, pcs));
            joints(n).pj = find(arrayfun(@(p) prod(joints(n).j-p.irange)<=0, pcs));
          otherwise
            ms = cell(joints(n).type,1);
            joints(n).ps = zeros(1,joints(n).type);
            for ii=1:joints(n).type
                ms{ii} = find(all(proccoords==origcoords(joints(n).is(ii),:), 2));
                if ii>1
                    ms{ii} = setdiff(ms{ii}, ms{ii-1});
                end
                joints(n).is(ii) = ms{ii}(1);
                joints(n).ps(ii) = find(arrayfun(@(p) prod(joints(n).is(ii)-p.irange)<=0, pcs));
            end
        end
    end
    
    % 4. Preprocess Excitation Information (and put into pcs)
    rlist = {'drcofsdw', 'drcofsdxi'};
    for n=1:length(excs)
        %         mi = find(all(proccoords==origcoords(excs(n).i,:),2));  % Doesn't
        %         work due to roundoff in certain versions
        mi = find(sum(abs(proccoords-origcoords(excs(n).i,:)),2)<1e-13);
        if length(mi)~=2
            error('Multiplicity of excitation point: %d. Unknown case.', length(mi));
        end
        if ~isfield(excs(n), 'drcofsdw') || isempty(excs(n).drcofsdw)
            excs(n).drcofsdw = matlabFunction(diff(excs(n).rcofs(w,xi), w), 'Vars', [w xi]);
        end
        if ~isfield(excs(n), 'drcofsdxi') || isempty(excs(n).drcofsdxi)
            excs(n).drcofsdxi = matlabFunction(diff(excs(n).rcofs(w,xi), xi), 'Vars', [w xi]);
        end
        % static forcing
        if ~isfield(excs(n), 'rcofs0') || isempty(excs(n).rcofs0)
            excs(n).rcofs0 = matlabFunction(excs(n).rcofs(0, xi), 'Vars', xi);  % Default behavior
        end
        if ~isfield(excs(n), 'drcofsdxi0') || isempty(excs(n).drcofsdxi0)
            excs(n).drcofsdxi0 = matlabFunction(diff(excs(n).rcofs0(xi), xi), 'Vars', xi);
        end

        if ~isfield(excs(n), 'nh') || isempty(excs(n).nh)
            excs(n).nh = 1;
        end

        excs(n).i = mi(1);
        excs(n).j = mi(2);

        pi = find(arrayfun(@(p) prod(excs(n).i-p.irange)<=0, pcs));  % Piece of excitation
        excs(n).pi = pi;

        pcs(pi).exci = [pcs(pi).exci; excs(n).i-pcs(pi).irange(1)+1];
        pcs(pi).excnh = [pcs(pi).excnh; excs(n).nh];
        pcs(pi).exccofs = {pcs(pi).exccofs{:}; 
                           {excs(n).rcofs, excs(n).drcofsdw, excs(n).drcofsdxi}};
        pcs(pi).exccofs0 = {pcs(pi).exccofs0{:}; 
                            {excs(n).rcofs0, excs(n).drcofsdxi0}};
    end
end
