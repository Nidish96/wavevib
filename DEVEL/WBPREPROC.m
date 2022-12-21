function [pcs, bcs, joints, excs] = WBPREPROC(pcs, bcs, joints, excs, wcomps)

    origcoords = cell2mat({pcs.coords}');
    if size(origcoords,2)<3
        origcoords = [origcoords zeros(size(origcoords,1),3)];
        origcoords = origcoords(:,1:3);
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

        pcs(i).N = size(pcs(i).coords,1);  % Number of Elements
        pcs(i).L = vecnorm(diff(pcs(i).coords([1 end],:)));
        pcs(i).nhat = diff(pcs(i).coords([1 end],:))/pcs(i).L;
        pcs(i).U = U;
        pcs(i).S = S;
        pcs(i).V = V;

        pcs(i).theta = atan2(vecnorm(pcs(i).nhat(1:2)), pcs(i).nhat(3));
        pcs(i).phi = atan2(pcs(i).nhat(2), pcs(i).nhat(1));
        pcs(i).irange = Nbef + [1 pcs(i).N];
        Nbef = Nbef + pcs(i).N;
    end
    proccoords = cell2mat({pcs.coords}');

    % 2. Preprocess Boundary Conditions
    Nwc = size(wcomps,1);
    dlist = {'dcofsdw', 'dcofsdxi'};
    for n=1:length(bcs)
        for di=dlist
            if ~isfield(bcs(n), di{:}) || isempty(bcs(n).(di{:}))
                bcs(n).(di{:}) = @(w,xi) zeros(1,Nwc);
            end
        end

        % Account for new nodes included in model
        mi = find(all(proccoords==origcoords(bcs(n).i,:), 2));
        bcs(n).i = mi(1);
        bcs(n).pi = find(arrayfun(@(p) prod(bcs(n).i-p.irange)<=0, pcs));
    end

    % 3. Preprocess Joints
    nlflist = {'dnlfcofsdw', 'dnlfcofsdxi'};
    nldlist = {'dnldcofsdw', 'dnldcofsdxi'};
    for n=1:length(joints)
        if ~isfield(joints(n), 'cofs') || isempty(joints(n).cofs)
            switch joints(n).type
                case {1, 2}
%                     joints(n).cofs = @(w,xi) [1 1 -1 -1; 1 -1 -1 1];  % Bar
                    joints(n).cofs = @(w,xi) [eye(Nwc) -eye(Nwc)];
                otherwise
                    error('Needs to be implemented still.');
            end
        end
        for di=dlist
            if ~isfield(joints(n), di{:}) || isempty(joints(n).(di{:}))
                joints(n).(di{:})= @(w,xi) zeros(Nwc,joints(n).type*Nwc);
            end
        end
        for di=nlflist
            if ~isfield(joints(n), di{:}) || isempty(joints(n).(di{:}))
                joints(n).(di{:})= @(w,xi) zeros(Nwc,1);
            end
        end
        for di=nldlist
            if ~isfield(joints(n), di{:}) || isempty(joints(n).(di{:}))
                joints(n).(di{:})= @(w,xi) zeros(1,joints(n).type*Nwc);
            end
        end
        if ~isfield(joints(n), 'nl')
            joints(n).nl = [];
        end

        % Account for new nodes included in model
        mi = find(all(proccoords==origcoords(joints(n).i,:), 2));
        mj = find(all(proccoords==origcoords(joints(n).j,:), 2));
        switch joints(n).type
            case {1, 2}
                joints(n).i = mi(1);
                mj = setdiff(mj, mi(1));
                joints(n).j = mj(1);
            otherwise
                error('Needs to be implemented still.');
        end
        joints(n).pi = find(arrayfun(@(p) prod(joints(n).i-p.irange)<=0, pcs));
        joints(n).pj = find(arrayfun(@(p) prod(joints(n).j-p.irange)<=0, pcs));
    end

    % 4. Preprocess Excitation Information (and put into pcs)
    rlist = {'drcofsdw', 'drcofsdxi'};
    for n=1:length(excs)
        mi = find(all(proccoords==origcoords(excs(n).i,:),2));
        if length(mi)~=2
            error('Multiplicity of excitation point: %d. Unknown case.', length(mi));
        end
        for di=rlist
            if ~isfield(excs(n), di{:}) || isempty(excs(n).(di{:}))
                excs(n).(di{:})= @(w,xi) zeros(Nwc,1);
            end
        end
        if ~isfield(excs(n), 'nh') || isempty(excs(n).nh)
            excs(n).nh = 1;
        end

        excs(n).i = mi(1);
        excs(n).j = mi(2);

        pi = find(arrayfun(@(p) prod(excs(n).i-p.irange)<=0, pcs));  % Piece of excitation
        excs(n).pi = pi;

        pcs(pi).exci = [pcs(pi).exci; excs(n).i];
        pcs(pi).excnh = [pcs(pi).excnh; excs(n).nh];
        pcs(pi).exccofs = {pcs(pi).exccofs{:}; 
            {excs(n).rcofs, excs(n).drcofsdw, excs(n).drcofsdxi}};
    end
end