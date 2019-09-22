% ECE310_Homework_1
% Kevin Kerliu
clear;
close all;
clc;
%%

% Question 4

% Part C
T = 2; % Period
N = 10; % Number of pulses

syms t
p(t) = triangularPulse(0,T/2,T,t);

x = 0;
for n=0:N-1
    x = x + p(t-n*T);
end
X = fourier(x);

% Graph x(t)
figure;
fplot (x, [0 N*T]);
title("x(t)");

% Graph |X(f)|
figure;
fplot(abs(X));
title("|X(f)|");

%%

% Question 6
% Sine wave amplitude
A = 2;
% Sine wave frequency
f = 20e3;
% Sine wave phase
theta = 0;

var = 0.2; % Variance
fs = 100e3; % Sampling rate
N = 1024; % N-pt DFT
n = 1000; % Samples

% Sine wave of amplitude 2 and frequency 20e3
t  = (0:n-1)/fs;
x = A*sin(2*pi*f*t + theta);
% AWGN of variance 0.2
noise = randn(1,n)*sqrt(var);
% Corrupted sine wave
x_corrupt = x + noise;
% Chebyshev window of length 100 and peak sidelobe level 30 dB
window = chebwin(1000, 30);
% Signal * Window Fucntion
x_windowed = x_corrupt.*window';
% Graph the magnitude spectrum
Yf = fftshift(fft(x_windowed,N));
Yf_dB = 20*log10(abs(Yf));
f_Hz = (-N/2:N/2-1)*(fs/N);
figure;
stem (f_Hz, Yf_dB);
xlabel("f (Hz)");
ylabel ("|Y(f)| (dB)");

%%

% Question 7
% x[n] = v[n] + 0.1v[n-1] - 0.72v[n-2] - 0.95x[n-1] -  0.9025x[n-2]
% H(z) = X(z)/V(z) 
%      = ( 1 + 0.1*z^-1 - 0.72*z^-2 ) / ( 1 + 0.95*z^-1 + 0.9025*z^-2 ) 
 
% Part A
% Numerator of H(z)
b = [1, 0.1, -0.72];
% Denominator of H(z)
a = [1, 0.95, 0.9025];
% Compute the poles and zeros of H(z)
[z,p,k] =  tf2zp(b,a);
% Plot the poles and zeros of H(z)
figure;
zplane(b,a);
title("Pole-Zero Plot of H(z)");

% Part B
% H(z) is the innovations filter of x
% x is ARMA.

% Part C
% Obtain an explicit formula for the power spectral density of x
% Sx(w) = |H(w)|^2 * Sv(w), where
% Sv(w) = 4 and
% H(w) = (1 + 0.1*e^-jw - 0.72*e^-jw2)/(1 + 0.95*e^-jw + 0.9025*e^-jw2)
% Sx(w) = 4 * |(1 + 0.1*e^-jw - 0.72*e^-jw2)/(1 + 0.95*e^-jw + 0.9025*e^-jw2)|^2

% Part D
% Generate N samples of v
N = 10^5;
% v[n] is mean 0, variance 4 white noise
% randn generates random numbers based on the Normal Distribution
% The default mean is 0 and the default variance is 1
v = sqrt(4)*randn(1, N);
% Apply the filter to v to generate N samples of x
x = filter(b, a, v);

% Part E
% Use the samples of x to estimate the PSD using the modified Welsh
% peridoogram
[s_est, w] = pwelch(x, hamming(512),256,512);
% Normalize the estimated PSD
s_est_normalized = s_est/mean(s_est);


% Part F
% H(w)
Hw = freqz(b,a,w);
% Sx(w)
Sx = 4*(abs(Hw)).^2;
% Normalize the exact PSD
Sx_normalized = Sx/mean(Sx);

% Part G
% Plot the PSDs
figure;
plot(w, Sx_normalized);
hold on; 
plot(w, s_est_normalized);
legend("Normalized Exact PSD", "Normalized Estimated PSD");
title("Exact PSD vs Estimated PSD");
xlabel("Normalized Digital Radian Frequency");
ylabel("PSD");
xlim([0, pi]);

% Part H
% Compute the angles of the poles of the transfer function
angle_peak_w = angle(p)
% Determine the maximum of the PSD
% Maximum of the estimated PSD
[peakPSD,maxWIndex] = max(s_est_normalized);
est_peak_w = w(maxWIndex)
% Maximum of the exact PSD
[peakPSD,maxWIndex] = max(Sx_normalized);
exact_peak_w = w(maxWIndex)
% The angles of the poles and the peak w are very close!