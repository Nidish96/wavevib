Nt = 32;

Om = 5;  % rad/s
% Setting Fundamental Frequency as (Om/3)
t = linspace(0, 2*pi/(Om/3), Nt+1)'; t = t(1:end-1);

yt = 45 + cos(Om*t) + 3*sin(Om*t) + cos(3*Om*t)*50;


% Setting the desired multiples of fundamental frequency ('harmonics')
h = [0; 1; 2; 3; 9];

% freq to time.
% format: [a0; a1; b1; a2; b2; ...]
Y = AFT(yt(:), h, Nt, 't2f');

% time to frequency
% format: a0+a1*cos(w*t)+b1*sin(w*t)+a2*cos(2*w*t)+b2*sin(2*w*t)+...
%   where 'w' is the fundamental frequency (set as Om/3 in line 5 above).
yt_recon = AFT(Y, h, Nt, 'f2t');