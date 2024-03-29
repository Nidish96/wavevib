function [yout] = AFT(yin, h, Nt, dir)
%AFT is the Generalized AFT routine that supports periodic as well as
%quasi-periodic inputs, supporting both time-frequency ('t2f') as well as
%frequency-time ('f2t') transformations.
%
%   USAGE :
%       youts = AFT(yin, h, Nt, dir);
%   INPUTS :
%       yin     : 't2f'(Nt^Nc, Ny) or 'f2t'(Nhc,Ny)
%       h       : (Nh,Nc)
%       Nt      :
%       dir     : 't2f' or 'f2t'
%   OUTPUT :
%       yout    : 't2f'(Nhc,Ny) or 'f2t'(Nt^Nc,Ny)

%%
    if sum(all(h==0, 2)) && ~all(h(1,:)==0)
        error('If you want the dc term, put it in the beginning of h')
    end

    Nc = size(h,2);
    Nh = size(h,1);
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
      
    iof = fix(Nt/2)+1;
    i0 = num2cell(repmat(iof, 1, Nc));
    
    Ny = size(yin,2);
    
    if strcmp(dir, 't2f')
        yout = zeros(Nhc, Ny);
        if Ny<Nc    % Option 1 - Loop over Ny
            for yi=1:Ny
                yf = reshape(yin(:, yi), [repmat(Nt, 1, Nc) ones(1, Nc==1)]);
                
                % Fourier Transform & Scale
                yf = fftshift(fftn(yf))*2/(Nt^Nc);
                yf(i0{:}) = yf(i0{:})/2;
        
                % Choose selected Harmonics
                k = 1;
                for hi=1:Nh
                    if all(h(hi,:)==0)  % zero harmonics
                        yout(k,yi) = yf(i0{:});
                        k = k+1;
                    else
                        id = num2cell(iof+h(hi,:)); % Focus for optimization
                        yout(k,yi)   = real(yf(id{:}));
                        yout(k+1,yi) = -imag(yf(id{:}));
                        k = k+2;
                    end
                end
            end
        else    % Option 2 - Do 1D fft over each direction
            yf = reshape(yin, [repmat(Nt, 1, Nc) Ny]);
            if Nc==2
                yf = fftshift(fftshift(fft2(yf),1),2);
            else
                for ic=1:Nc
                    yf = fftshift(fft(yf, [], ic),ic);
                end
            end
            yf = yf*2/(Nt^Nc);
            yf(i0{:}, :) = yf(i0{:}, :)/2;
            % Choose Selected Harmonics
            k = 1;
            for hi=1:Nh
                if all(h(hi,:)==0)
                    yout(k,:) = yf(i0{:}, :);
                    k = k+1;
                else
                    id = num2cell(iof+h(hi,:));
                    yout(k, :) = real(yf(id{:}, :));
                    yout(k+1, :) = -imag(yf(id{:}, :));
                    k = k+2;
                end
            end
        end
    elseif strcmp(dir, 'f2t')
        yout = zeros(Nt^Nc,Ny);
        yf = zeros([repmat(Nt, 1, Nc) Ny]);
        k = 1;
        for hi=1:Nh
            if all(h(hi,:)==0)
                yf(i0{:},:) = yin(k,:)*2;
                k = k+1;
            else
                id = num2cell(iof+h(hi,:));
                yf(id{:}, :) = yin(k,:)-1j*yin(k+1,:);
                id = num2cell(iof-h(hi,:));
                yf(id{:}, :) = yin(k,:)+1j*yin(k+1,:);
                k = k+2;
            end
        end
        yf = yf*(Nt^Nc)/2;

        if Ny<Nc
            % IFFT
            for yi=1:Ny  % can be task-parallelized
                yout(:, yi) = reshape(ifftn(ifftshift(yf)), Nt^Nc, 1);
            end
        else    % Option 2 - Do 1D fft over each direction
            if Nc==2
                yf = ifft2(ifftshift(ifftshift(yf,1),2));
            else
                for ic=1:Nc
                    yf = ifft(ifftshift(yf,ic), [], ic);
                end
            end
            yout = real(reshape(yf, Nt^Nc, Ny));
        end
    end
end