function fig2gif(filename, LoopCount, varargin )
%% DANIEL FOCHLER ILA 2016
%% BESCHREIBUNG
%   Die Funktion dient dazu in Matlab eine Figure, welche Schritt für Schritt
%   aufgebaut wird.
%   Anwendung siehe Beispiel:
%   Beachte die Befehle: drawnow (notwendig); cla (optional)
%% INPUT
%  siehe Beispiel
%% BEISPIEL
% x = 0:0.01:1;
% filename = 'testnew.gif';
% a = 1:0.5:5;
% for jj = 1:length(a);
%       figure('units','normalized','outerposition',[0 0 1 1], 'Color', white(1));
%       y = a(jj).*x;
%       plot(x,y)
%       drawnow
%       fig2gif(filename, jj, .2)
%       close gcf
% end

%% Inputs
if nargin >= 3
    dt = varargin{1};
else
    dt = .1;
end
if nargin >= 4
    hand = varargin{2}; 
else 
    hand  = gcf;
end

%%
[fpth,fnm,fex]=fileparts(filename);
if ~strcmp(fex(2:end), '.gif')
    filename = [fpth, '/', fnm, '.gif'];
end
%%
frame = getframe(hand);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if LoopCount == 1
    imwrite(imind,cm,filename,'gif', 'LoopCount',inf, 'DelayTime', dt);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', dt);
end
%%
