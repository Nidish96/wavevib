function [uout] = AFT(uin, Nt, h, dir)
%AFT conducts time-to-frequency & frequency-to-time transforms
%
%   USAGE :
%       uout = AFT(uin, Nt, h, dir);
%   INPUTS :
%       uin     : (Nt, Np) time-series of "Np" signals
%                       (OR)
%                 (Nhc, Np) harmonic coefficients of "Np" signals
%       Nt      : Number of time points
%       h       : List of harmonic coefficients
%       dir     : 't2f' time to frequency transform
%                       (OR) 
%                 'f2t' frequency to time transform
%   OUTPUTS :
%       uout    : (Nhc, Np) harmonic coefficients of the "Np" signals
%                       (OR)
%                 (Nt, Np) time-series of the "Np" signals

    Nhc = sum((h==0)+2*(h~=0));
    Np = size(uin, 2);
    if strcmp(dir, 't2f')  % time to frequency transform
        uout = zeros(Nhc, Np);
        
        uh = fft(uin)/(Nt/2); 
        uh(1,:) = uh(1,:)/2;
        
        if h(1)==0  % zero harmonic is present in list
            uout(1, :) = uh(1, :);  % zero harmonic
            uout(2:2:end, :) = real(uh(1+h(2:end), :));  % cosine coefficients
            uout(3:2:end, :) = -imag(uh(1+h(2:end), :));  % sine coefficients
        else
            uout(1:2:end, :) = real(uh(1+h, :));  % cosine coefficients
            uout(2:2:end, :) = -imag(uh(1+h, :)); % sine coefficients
        end
    else  % frequency to time transform        
        uh = zeros(Nt, Np);
        if h(1)==0  % zero harmonic is present in list
            uh(1, :) = uin(1, :)*2;
            uh(1+h(2:end), :) = uin(2:2:end,:)-uin(3:2:end,:)*1j;
            uh(Nt-h(2:end)+1, :) = uin(2:2:end,:)+uin(3:2:end,:)*1j;
        else
            uh(1+h, :) = uin(1:2:end,:)-uin(2:2:end,:)*1j;
            uh(Nt-h+1, :) = uin(1:2:end,:)+uin(2:2:end,:)*1j;
        end
        uout = ifft(uh*(Nt/2));
    end
end
