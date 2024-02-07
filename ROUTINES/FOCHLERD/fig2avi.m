function fig2avi(filename, TimeStepCounter, NumberOfFrames,  varargin )
%% DANIEL FOCHLER ILA 2016
%% BESCHREIBUNG
% Die Funktion dient dazu in Matlab eine Figure, welche Schritt für Schritt
% aufgebaut wird.
% Anwendung siehe Beispiel:
% Beachte die Befehle: drawnow (notwendig); 
%% INPUT
%  siehe Beispiel
%% BEISPIEL
% x = 0:0.01:1;
% filename = 'testnew.avi';
% a = 1:0.5:5;
% for jj = 1:length(a);
%       figure('units','normalized','outerposition',[0 0 1 1], 'Color', white(1));
%       y = a(jj).*x;
%       plot(x,y)
%       drawnow
%       fig2avi(filename,jj, length(a), .2)
%       close gcf
% end

%% Inputs
    if nargin >= 4
        dt = varargin{1};
    else
        dt = .1;
    end
    if nargin >= 5
        Quality = varargin{2};
    else
        Quality = 100;
    end

    %%
    frame = getframe(gcf);
    save(['./frame', num2str(TimeStepCounter)], 'frame');

    if TimeStepCounter == NumberOfFrames
        %% Make avi file
        [fpth,fnm,fex]=fileparts(filename);
        if ~strcmp(fex(2:end), '.avi')
            % Delete a previously generated file:
            moviefile = fullfile(filename);
            if exist(moviefile, 'file')
                delete(moviefile);
            end
        end
        % Create avi frame by frame
        
        writerObj = VideoWriter(moviefile);
        writerObj.FrameRate = 1/dt;
        writerObj.Quality = Quality;
        open(writerObj);
        
        for jj = 1:NumberOfFrames
            load(['frame', num2str(jj)]);
            writeVideo(writerObj, frame);
            delete(['frame', num2str(jj) '.mat']);
        end
        close(writerObj);
    end
end
