% ECE310_Homework_6
% Kevin Kerliu
close all;
clear;
clc;
%% Question 4

% Part A
% Define H0, H1, F0, F1

N = 12;

% H0
H0 = zeros(1,N);
H0(1) = 0.15774243;
H0(2) = 0.69950381;
H0(3) = 1.06226376;
H0(4) = 0.44583132;
H0(5) = -0.31998660;
H0(6) = -0.18351806;
H0(7) = 0.13788809;
H0(8) = 0.03892321;
H0(9) = -0.04466375;
H0(10) = -7.83251152e-4;
H0(11) = 6.75606236e-3;
H0(12) = -1.52353381e-3;

% H1
H1 = zeros(1,N);
for k=1:N
H1(k) = (-1)^(k+1) * H0(N+1-k);
end

% F0
F0 = zeros(1,N);
for k=1:N
F0(k) = H0(N+1-k);
end

% F1
F1 = zeros(1,N);
for k=1:N
F1(k) = -(-1)^(k+1) * H0(k);
end

%%

% Part B

% Define T and A
T = .5 * (conv(F0,H0) + conv(F1,H1));
A = .5 * (conv(F0,-F1) + conv(F1,F0)); 
% f1[n] = -(-1)^n * h0[n] 
% f0[n] = (-1)^n * h1[n]
% so we can make the substitutions

% c = 2, N = 11
% We find c at the 12th index of T
% This corresponds to an offset of 11

% Determining error (deviation from expected value of zero)
errA = max(abs(A)); 
nonzeroTindex = find(abs(T) < max(abs(T)));
errT = max(abs(T(nonzeroTindex))); 
% errA = 5.551115123125783e-17
% errT = 0.001095774507172

% Part C
[freqzH0,wH0] = freqz(H0,1,1024);
[freqzH1,wH1] = freqz(H1,1,1024);
linMagH1 = abs(freqzH0);
linMagH2 = abs(freqzH1);
figure;
plot(wH0,linMagH1);
hold on;
plot(wH1,linMagH2);
xlim([0 pi]);

% Labeling
legend("H1","H2");
title("Magnitude Response of H0 and H1 (Linear Scale)");
xlabel("Frequency (rad/sec)");
ylabel("|H(w)| (Linear Scale)");

%%

% Part D

% Determining the pcc and error (deviation from pcc)
H0H1sumsquares = abs(freqzH0).^2 + abs(freqzH1).^2;
pcc = round(H0H1sumsquares); % power complementary constant
pcc = pcc(1);
maxError = max(pcc - min(H0H1sumsquares),pcc - max(H0H1sumsquares));

% pcc = 4
% maxError = 0.0063
%%

% Part E
% Define E = [E00, E01; E10, E11]
E00 = H0(1:2:end);
E01 = H0(2:2:end);
E10 = H1(1:2:end);
E11 = H1(2:2:end);

% Define P = paraconj(E) * E
P00 = conv(flip(E00),E00) + conv(flip(E10),E10);
P01 = conv(flip(E00),E01) + conv(flip(E10),E11);
P10 = conv(flip(E01),E00) + conv(flip(E11),E10);
P11 = conv(flip(E01),E01) + conv(flip(E11),E11);

% Determining error (deviation from expected value of zero)
errP01 = max(abs(P01));
errP10 = max(abs(P10));
nonzeroP00index = find(abs(P00) < max(abs(P00)));
errP00 = max(abs(P00(nonzeroP00index)));
nonzeroP11index = find(abs(P11) < max(abs(P11)));
errP11 = max(abs(P11(nonzeroP11index)));
d_diff = abs(max(abs(P00)) - max(abs(P11)));

% errP01 = 2.775557561562891e-17
% errP10 = 2.775557561562891e-17
% errP00 = 0.001095774507172
% errP11 = 0.001095774507172
% d_diff = 2.220446049250313e-16
% d = 2 at n0 = 0
