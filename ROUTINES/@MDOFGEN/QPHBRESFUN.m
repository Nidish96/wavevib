function [R, dRdU, FNL] = QPHBRESFUN(m, Up, ws, Fl, h, Nt, tol, varargin)
%QPHBRESFUN 
%
%   USAGE: 
%       [R, dRdU, dRdw, FNL] = QPHBRESFUN(m, Up, ws, Fl, h, Nt, tol)
%   INPUTS:
%       MDOFGEN class
%       Up
%       Fl
%       h
%       Nt
%       tol 
%   OUTPUTS:
%       R
%       dRdU
%       dRdw
%       FNL
    
    Nc = size(h,2);  % Number of components
    Nh = size(h,1);  % Number of harmonics
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));  % Number of harmonic coefficients
    
    Ws = Up(end)*ws;  % real frequency is ws scaled by Up(end)

    E = QPHARMONICSTIFFNESS(m.M, m.C, m.K, Ws, h);  % Harmonic Stiffness
    
    tau = linspace(0, 2*pi, Nt+1); tau(end) = [];
    taus = cell(Nc, 1);
    [taus{:}] = ndgrid(tau);
    t = zeros(repmat(Nt, 1, Nc));
    for ic=1:Nc
        t = t+taus{ic}*Ws(ic);
    end
    
    D1 = QPHARMONICSTIFFNESS(0, 1, 0, Ws, h);  % time derivative matrix
    
    cst = QPTIMETRANS(eye(Nhc), h, Nt);  % basis functions
    sct = QPTIMETRANS(D1, h, Nt);  % basis function derivatives
    
    FNL = zeros(m.Ndofs*Nhc, 1);
    dFNL = zeros(m.Ndofs*Nhc);
    for ni=1:length(m.NLTs)
        Unl = (m.NLTs(ni).L*reshape(Up(1:end-1), m.Ndofs, Nhc))';  % Nhc x Ndnl
        Ndnl = size(m.NLTs(ni).L, 1);

        unlt   = QPTIMETRANS(Unl, h, Nt);
        unldot = QPTIMETRANS(D1*Unl, h, Nt);
        
        if mod(m.NLTs(ni).type, 2)==0  % INSTANTANEOUS FORCING
            [ft, dfdu, dfdud] = m.NLTs(ni).func(taus, unlt, unldot);
        
            F = QPFOURIERCOEFF(ft, h);
            J = zeros(size(m.NLTs(ni).L,1)*Nhc, size(m.NLTs(ni).L,1)*Nhc);
            dFdU = reshape(QPFOURIERCOEFF(reshape(dfdu.*permute(cst, [1, 3, 2]), Nt^Nc, Ndnl*Nhc), h) + ...
                QPFOURIERCOEFF(reshape(dfdud.*permute(sct, [1, 3, 2]), Nt^Nc, Ndnl*Nhc), h), ...
                Nhc, Ndnl, Nhc);
        
            for di=1:Ndnl
                J(di:Ndnl:end, di:Ndnl:end) = dFdU(:, di, :);
            end
        else % HYSTERETIC FORCING
            Nnlf = size(m.NLTs(ni).L, 1);
            if ~isempty(m.NLTs(ni).Lf)
                Nnlf = size(m.NLTs(ni).Lf, 2);
            end
            if Nnlf>1 || size(unlt, 2)>1
                error('Not implemented for multi-DOF dependent multi-force nonlinearities yet')
            end

            switch m.NLTs(ni).qptype
                case {1,2}  % Option 1: Solve using NSOLVE or fsolve
                    % Construct Nmat
                    Nmat = CONSTRUCTNMAT(Ws, Nc, Nt, m.NLTs(ni).qptype);
                    ft = ones(Nt^Nc, Nnlf);
                    
                    opt = struct('Display', false);
        %             tic
                    ft = NSOLVE(@(ft) QPMARCHRESFUN(ft, unlt, ...
                        @(u, up, fp) m.NLTs(ni).func(0, u, 0, 0, up, fp), ...
                        Nmat), ft, opt);
        %             toc
        %             opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'off');
        %             ft = fsolve(@(ft) QPMARCHRESFUN(ft, unlt, ...
        %                 @(u, up, fp) m.NLTs(ni).func(0, u, 0, 0, up, fp), ...
        %                 Nmat), ft, opt);
        
                    % Get Jacobians
                    [~, dresdf, dresdu] = QPMARCHRESFUN(ft, unlt, ...
                        @(u, up, fp) m.NLTs(ni).func(0, u, 0, 0, up, fp), ...
                        Nmat);
                    dfdut = -dresdf\dresdu;
                    
                    F = QPFOURIERCOEFF(ft, h);
                    J = QPFOURIERCOEFF(dfdut*cst, h);
                case 3  % Option 2: Sequential Solution: Only possible for Nmtype=3
                    [Nmat, bpis, bpjs] = CONSTRUCTNMAT(Ws, Nc, Nt, m.NLTs(ni).qptype);
                    bpjs = cellfun(@(c) unique(c), bpjs, 'UniformOutput', false);
        
                    its = 0;
                    ft = zeros(Nt^Nc, Nnlf);
                    dfdai = zeros(Nt^Nc, Nhc);
                    fprev = ft(bpis{1},:);
        %             tic
                    while its==0 || max(abs(ft(bpis{1})-fprev))>tol 
                        fprev = ft(bpis{1},:);
                        for ti=1:Nt
                            % Only force estimation
                            ft(bpis{ti}) = m.NLTs(ni).func(0, unlt(bpis{ti}), 0, 0, Nmat(bpis{ti},bpjs{ti})*unlt(bpjs{ti}), Nmat(bpis{ti},bpjs{ti})*ft(bpjs{ti}));
        
        %                     % Force & Jacobian estimation - Naive version
        %                     [ft(bpis{ti}), ~, dfdfp, dfdu, dfdup] = m.NLTs(ni).func(0, unlt(bpis{ti}), 0, 0, Nmat(bpis{ti},:)*unlt, Nmat(bpis{ti},:)*ft);
        %                     dfdai(bpis{ti},:) = diag(dfdfp)*Nmat(bpis{ti},:)*dfdai + dfdu.*cst(bpis{ti},:) + diag(dfdup)*Nmat(bpis{ti},:)*cst;
        
        %                     % Force & Jacobian estimation - respectful of sparsity in Nmat
        %                     [ft(bpis{ti}), ~, dfdfp, dfdu, dfdup] = m.NLTs(ni).func(0, unlt(bpis{ti}), 0, 0, Nmat(bpis{ti},bpjs{ti})*unlt(bpjs{ti}), Nmat(bpis{ti},bpjs{ti})*ft(bpjs{ti}));
        %                     dfdai(bpis{ti},:) = diag(dfdfp)*Nmat(bpis{ti},bpjs{ti})*dfdai(bpjs{ti},:) + dfdu.*cst(bpis{ti},:) + diag(dfdup)*Nmat(bpis{ti},bpjs{ti})*cst(bpjs{ti},:);
                        end
                        its = its+1;
                    end
                    % Iterate and just estimate the Jacobians in the end
                    for ti=1:Nt  % Run the iterations once more
                        [ft(bpis{ti}), ~, dfdfp, dfdu, dfdup] = m.NLTs(ni).func(0, unlt(bpis{ti}), 0, 0, Nmat(bpis{ti},bpjs{ti})*unlt(bpjs{ti}), Nmat(bpis{ti},bpjs{ti})*ft(bpjs{ti}));
                        dfdai(bpis{ti},:) = diag(dfdfp)*Nmat(bpis{ti},bpjs{ti})*dfdai(bpjs{ti},:) + dfdu.*cst(bpis{ti},:) + diag(dfdup)*Nmat(bpis{ti},bpjs{ti})*cst(bpjs{ti},:);
                    end
%                     fprintf('%d\n',its)
                    F = QPFOURIERCOEFF(ft, h);
                    J = QPFOURIERCOEFF(dfdai, h);
        %             toc
            end
            
%             J = zeros(size(m.NLTs(ni).L,1)*Nhc, size(m.NLTs(ni).L,1)*Nhc);
%             for di=1:Ndnl
%                 for dj=1:Ndnl
%                     tmp = squeeze(dfdu(:, di, dj, :));
%                     if ~isempty(find(tmp~=0, 1))
%                         J(di:Ndnl:end, dj:Ndnl:end) = ...
%                             GETFOURIERCOEFF(h, tmp);
%                     end
%                 end
%             end
        end
        
        if m.NLTs(ni).type<=5  % Self adjoint forcing
            FNL  = FNL + reshape(m.NLTs(ni).L'*F', Nhc*m.Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).L')*J*kron(eye(Nhc), m.NLTs(ni).L);
        else  % Non-self adjoint forcing
            FNL  = FNL + reshape(m.NLTs(ni).Lf*F', Nhc*m.Ndofs, 1);
            dFNL = dFNL + kron(eye(Nhc), m.NLTs(ni).Lf)*J*kron(eye(Nhc), m.NLTs(ni).L);                
        end
    end
    
    % Residue
    if length(m.Rsc)~=length(Fl)
        m.Rsc = (1/max(abs(Fl)))*ones(length(Fl),1);
    end
    R = [E*Up(1:end-1) + FNL - Fl].*m.Rsc;
    dRdU = (E+dFNL).*m.Rsc;
%     dRdw = (dEdw*Up(1:end-1)).*m.Rsc;

    % All Gradient terms in one matrix
%     if ~isempty(varargin)
%         dRdU = [dRdU dRdw];
%     end
end

%%
function [res, dresdf, dresdu] = QPMARCHRESFUN(ft, unlt, func, Nmat)
%Residue function for 
    [fnl, ~, dfnldfp, dfnlduc, dfnldup] = func(unlt, Nmat*unlt, Nmat*ft);
    
    res = ft-fnl;
    dresdf = sparse(eye(length(ft))-diag(dfnldfp)*Nmat);
    dresdu = sparse(-diag(dfnlduc) - diag(dfnldup)*Nmat);
%             [ft(ti,:), dfdu(ti,:,:,:)] = ...
%                 m.NLTs(ni).func(t(ti), unlt(ti,:), h, t(tm1), ...
%                 unlt(tm1,:), ft(tm1,:), dfdu(tm1,:,:,:));    
end
%%
function [Nmat, varargout] = CONSTRUCTNMAT(ws, Nc, Nt, varargin)  % Need to think about w-gradients
    if nargin==3
        Nmtype = 1;
    else
        Nmtype = varargin{1};
    end
    deltau = 2*pi/Nt; 
    taus = cell(Nc, 1);
    [taus{:}] = ndgrid(1:Nt);    
    switch Nmtype
        case 1
            dt_vec = ws*deltau/vecnorm(ws);  % vector corresponding to deltatau amplitude in real time dxn

            ptsb = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^Nc-1))));  % using binary for the construction of points on a unit square
            oppi = (2^Nc):-1:1;    % diagonally opposite points are retrieved using the binary inverses
            xis = ptsb*(-deltau);  % coordinates in "tau" space relative to origin

            Lm = deltau^Nc;                                 % Lebesgue Measure of each cell in tau space
            Nsf = prod(abs(xis(oppi,:)-(-dt_vec)), 2)'/Lm;  % Lagrange Shape Functions to interpolate previous point in the cell diagonally behind (in tau space)

            ijs = fix(cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)'));  % indices of all points
            evijs = mod(repmat(ijs, 1, 1, 2^Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;   % indices of points forming the cell that is diagonally behind
            evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);         % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind

            % Build sparse interpolation matrix
            Nmat = sparse(repmat((1:Nt^Nc)', 1, 2^Nc), evns, repmat(Nsf, Nt^Nc, 1));
        case 2
            ptsb = eye(Nc);

            Nsf = ws/sum(ws);  % Shape functions constructed using omegas

            ijs = cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)');  % indices of all points
            evijs = mod(repmat(ijs, 1, 1, Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the FD stencil
            evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points from the stencil

            % Build sparse interpolation matrix
            Nmat = sparse(repmat((1:Nt^Nc)', 1, Nc), evns, repmat(Nsf, Nt^Nc, 1));
        case 3
            ptsb = ones(1, Nc);  % diagonally opposite point
            
            ijs = cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)');  % indices of all points
            evijs = mod(repmat(ijs, 1, 1, 1) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the cell that is diagonally behind
            
            evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind
            
            % Build interpolant matrix
            dt_vecs = ws * deltau./ws(:);  % Each row is the "previous point" projected on a different plane. We will have to apply the correct one for the correct points
            
            ptsbm = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^(Nc-1)-1))));  % using binary for the construction of points on a unit square
            Lm = deltau^(Nc-1);
            
            inds = reshape((1:Nt^Nc)', [repmat(Nt, 1, Nc) ones(Nc==1)]);
            S.subs = repmat({':'}, 1, Nc);
            S.type = '()';
            bpis = cell(Nt,1);  bpjs = cell(Nt,1);    vals = cell(Nt,1);
            for ti=1:Nt  % march over the diagonal
                bpis{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,1);
                bpjs{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,2^(Nc-1));
                vals{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,2^(Nc-1));
                for ci=1:Nc
                    S.subs = repmat({ti:Nt}, 1, Nc);
                    S.subs{ci} = ti;
                    sinds = unique(reshape(subsref(inds, S), [], 1));
                
                    bpis{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1))) = sinds;
            
                    % build local Lagrange interpolant for these points
                    cinds = setxor(1:Nc, ci);
                    dt_loc = -dt_vecs(ci, cinds);
            
                    n1s_loc = floor(dt_loc/deltau);  % relevant point on the diagonal
                    
                    xis = (n1s_loc+ptsbm)*deltau;  % points on a cube around previous projection
            
                    Nsfs = abs(prod(xis(end:-1:1,:)-dt_loc,2)'/Lm);  % Lagrangian Shape Functions
                    
                    ptsb = zeros(2^(Nc-1), Nc);
                    ptsb(:, cinds) = n1s_loc+ptsbm+1;
            
                    bpevijs = mod(repmat(ijs(evns(sinds),:), 1, 1, 2^(Nc-1)) + permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the relevant cell
                    bpjs{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1)),:) = squeeze(sum((Nt.^(0:Nc-1)).*(bpevijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) 
                    vals{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1)),:) = repmat(Nsfs, (Nt-ti+1)^(Nc-1), 1);
                end
                [bpis{ti}, indtis] = unique(bpis{ti});
                bpjs{ti} = bpjs{ti}(indtis,:);
                vals{ti} = vals{ti}(indtis,:);
            end
            varargout{1} = bpis;
            varargout{2} = bpjs;
            [bpis, uinds] = unique(cell2mat(bpis));
            bpjs = cell2mat(bpjs);  bpjs = bpjs(uinds,:);
            vals = cell2mat(vals);  vals = vals(uinds,:);
            
            % Build sparse interpolation matrix
            Nmat = sparse(repmat(bpis, 1, 2^(Nc-1)), bpjs, vals);
    end
end