function [] = DEPICTBEAM_TM3D(Les,wd1s,wd2s,XYZc,U, varargin)
    
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

    plot3(XYZc(:,1), XYZc(:,2), XYZc(:,3), 'k.--'); hold on
    plot3(XYZc(:,1)+U(1:6:end), XYZc(:,2)+U(2:6:end), XYZc(:,3)+U(3:6:end), 'o-')
    for e=1:length(Les)
        xis = (e-1)*6 + [1 7];

        XYZce = XYZc(e+(0:1),:);
        XYZce(:,2) = XYZce(os,2);  XYZce(:,3) = XYZce(os,3);
%         DRAWCUBOID([Les(e); wd1s(e); wd2s(e)], ...
%            (XYZce+[U(xis) U(xis+1) U(xis+2)])', ...
%            [atan(U(xis+3)) atan(U(xis+4)) U(xis+5)]', col, alph);
        DRAWCUBOID([Les(e); wd1s(e); wd2s(e)], ...
           (XYZce+[U(xis) U(xis+1) U(xis+2)])', ...
           [atan(-U(xis+5)) atan(-U(xis+4)) atan(-U(xis+3))]', col, alph);
    end
end
