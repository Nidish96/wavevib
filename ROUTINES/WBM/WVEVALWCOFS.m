function [Us] = WVEVALWCOFS(aris, w, dL, Xs, pcs, Klib)
%WVEVALWCOFS This is a routine for evaluating the wave coefficients at
%arbitrary points on the model. 
%
%       CURRENTLY UNDER DEVELOPMENT

    if size(Xs,2)>1
        error('Need to implement 3D space coordinates')
    end
    Us = zeros(size(dL,1), length(Xs));
    Nwc = size(pcs(1).wcomps,1);
    for n=1:length(pcs)
        Ks   = pcs(n).wcomps(:,1).*cellfun(@(c) c(w,0), {Klib(pcs(n).wcomps(:,2)).K}.');

        for i=1:size(pcs(n).coords,1)-1
            inds = find(pcs(n).coords(i,1)<=Xs & pcs(n).coords(i+1,1)>=Xs);
            
            dX = Xs(inds)-pcs(n).coords(i,1);
            k = pcs(n).irange(1)+i-1;
            Us(:,inds) = dL*(aris((k-1)*Nwc+(1:Nwc)).*exp(Ks(:).*dX(:)'));
        end
    end
end