function [ F, Pxx ] = myPSD(x,Fsampling,N)
% Power Spectral Density
% [ F, Pxx ] = myPSD(x, Finterest, Fsampling)

x = x-mean(x);
% Finterest (Hz), the lowest freq we want to resolve accurately

% WindowSize = 32 * Fsampling / Finterest;
% NumOfWindows = 128;%length(x)/WindowSize;

WindowLength = N;%length(x) / NumOfWindows;
overlap = WindowLength / 2;

[pxx , f] = pwelch(x, WindowLength, overlap, Fsampling/2, Fsampling);

Pxx = pxx;
F = f;

end
