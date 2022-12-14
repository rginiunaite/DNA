%% fft
% run explore before hand

%x=(1:147);
%umaxall = sin(x);

Fs = 147;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 147;             % Length of signal


Y = fft(mean(umaxall));

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);


figure
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")

% period is Fs/(x value at the peak)


