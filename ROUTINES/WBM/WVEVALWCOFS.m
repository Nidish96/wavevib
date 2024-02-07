function [Us, Xs] = WVEVALWCOFS(ac, w, h, dL, dL0, Nx, pcs, Klib, varargin)
%WVEVALWCOFS This is a routine for evaluating the wave coefficients at
%arbitrary points on the model. 
%
%   USAGE:
%       [Us, Xs] = WVEVALWCOFS(ac, w, h, dL, dL0, Nx, pcs, Klib);
%   INPUTS:
%       ac       : (Npts*Nwc, Nh) vector of complex wave coefficients
%       w        : (1,1) frequency
%       h        : (Nh,1) list of harmonics
%       dL       : (1,Nwc) component summing vector
%       dL0      : (1,Nwc) static component summing vector    
%       Nx       : (1,1) Number of points per piece
%       pcs      : Pieces structure
%       Klib     : K library structure
%   OUTPUTS:
%       Us       : (Npcs, 1) cell, each with (Nx, Nh) wave coefficients
%       Xs       : (Npcs, 1) cell, each with (Nx, 3) coordinate points    

    if isempty(varargin)
        dero = 0;
    else
        dero = varargin{1};
    end
    
    Nh = length(h);

    Xs = arrayfun(@(pc) cell2mat(arrayfun(@(i) linspace(pc.coords(1,i), pc.coords(end,i), Nx)', 1:3, 'UniformOutput', false)), pcs, 'UniformOutput', false);
    Us = arrayfun(@(pc) zeros(Nx,Nh), pcs, 'UniformOutput', false);
    
    for hi=1:Nh
        for n=1:length(pcs)
            Ks = pcs(n).wcomps(:,1).*...
                 cellfun(@(c) c(h(hi)*w,0), {Klib(pcs(n).wcomps(:,2)).K}.');
            Nwc = length(Ks);
            
            x1 = (Xs{n}-pcs(n).coords(1,:))*pcs(n).V/pcs(n).S;  % Ordered 1D version of pts
            for i=1:pcs(n).N-1
                inds = find(pcs(n).U(i)-x1<=eps & x1-pcs(n).U(i+1)<=eps);
                dx = (x1(inds)-pcs(n).U(i))*pcs(n).S;
                
                k = pcs(n).irange(1)-1+i;
                if all(h(hi,:)~=0)
                    Us{n}(inds,hi) = dL*(ac((k-1)*Nwc+(1:Nwc),hi).*...
                                         ((Ks.^dero).*exp(Ks.*dx')) );
                else
                    % switch dero
                    %   case 0
                    %     Us{n}(inds,hi) = dL0*(ac((k-1)*Nwc+(1:Nwc),hi).*(dx.^(0:Nwc-1))');
                    %   otherwise
                        tmp = ac((k-1)*Nwc+(1:Nwc),hi).*...
                              (gamma(1:Nwc)./gamma((1:Nwc)-dero).*...
                               dx.^((0:Nwc-1)-dero))';
                        tmp(1:min(dero,Nwc), :) = 0;
                        Us{n}(inds, hi) = dL0*tmp;
                    % end
                end
            end
        end
    end
end
