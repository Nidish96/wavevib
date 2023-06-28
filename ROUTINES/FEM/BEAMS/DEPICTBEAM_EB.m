function [] = DEPICTBEAM_EB(Les,wd1s,wd2s,Xs,XYZc,U, varargin)
    
    col  = [0.5 0.5 0.5];
    alph = 0.5;
    os   = 1;
	if nargin>=7
        col = varargin{1};
    end
    if nargin>=8
        alph = varargin{2};
    end
    if nargin>=9
        os = varargin{3};
    end

    plot3(Xs, Xs*0, Xs*0, 'k.--'); hold on
    plot3(Xs+U(1:5:end), U(2:5:end), U(4:5:end), 'o-')
    for e=1:length(Les)
        xis = (e-1)*5 + [1 6];

        XYZce = XYZc(e+(0:1),:);
        XYZce(:,2) = XYZce(os,2);  XYZce(:,3) = XYZce(os,3);
        DRAWCUBOID([Les(e); wd1s(e); wd2s(e)], ...
           (XYZce+[U(xis) U(xis+1) U(xis+3)])', ...
           [atan(-U(xis+2)) atan(U(xis+4)) U(xis+4)*0]', col, alph);
    end
end