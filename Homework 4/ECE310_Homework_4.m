% ECE310_Homework_4
% Kevin Kerliu
clear;
close all;
clc;
%%

% Question 4 
[z,p,k] = ellip(3,2,30,[0.2,0.4]);
[b,a] = zp2tf(z,p,k);
[d0,d1] = tf2ca(b,a);
p0 = fliplr(d0);
p1 = fliplr(d1);

% Part A
[b_computed,a_computed] = compute(p0,d0,p1,d1);

% Part B
[H1,W1] = freqz(b,a);
[H2,W2] = freqz(b_computed, a_computed);
max_abs_error = max(abs(H2-H1));

% Part C
figure;

subplot(2,1,1);
magH1 = 20*log10(abs(H1));
plot(W1/pi,magH1);

title("Magnitude Response");
xlabel("Normalized Frequency (\times \pi rad/sample)");
ylabel("Magnitude (dB)");

ax = gca;
ax.YLim = [-50, 10];
ax.XTick = 0:.1:1;

subplot(2,1,2);
phasez(b,a);

% Part D
bscale = 8;

a_nearest_16 = round((a)*16)/(16);
b_nearest_16 = round((b*bscale)*16)/(bscale*16);

a_nearest_4 = round((a)*4)/(4);
b_nearest_4 = round((b*bscale)*4)/(bscale*4);

[h_nearest_16, w_nearest_16] = freqz(b_nearest_16,a_nearest_16);
[h_nearest_4, w_nearest_4] = freqz(b_nearest_4,a_nearest_4);

figure;

plot(W1/pi,20*log10(abs(H1)));

hold on;
plot(w_nearest_16/pi,20*log10(abs(h_nearest_16)));

hold on;
plot(w_nearest_4/pi,20*log10(abs(h_nearest_4)));

title("Magnitude Response");
xlabel("Normalized Frequency (\times \pi rad/sample)");
ylabel("Magnitude (dB)");
legend("Original", "Nearest 16th", "Nearest 4th");

ax = gca;
ax.YLim = [-50, 10];
ax.XTick = 0:.1:1;

% Part E
figure;
sgtitle("Pole-Zero Plots");

subplot(3,1,1);
zplane(b,a);
title("Original");

subplot(3,1,2);
zplane(b_nearest_16,a_nearest_16);
title("Nearest 16th");

subplot(3,1,3);
zplane(b_nearest_4,a_nearest_4);
title("Nearest 4th");

% Nearest 16th: Poles get closer to the unit circle.
% Nearest 4th: Some poles are outside the unit circle.

%  Part F
d0_nearest_4 = round((d0)*4)/(4);
d1_nearest_4 = round((d1)*4)/(4);
p0_nearest_4 = round((p0)*4)/(4);
p1_nearest_4 = round((p1)*4)/(4);
[b_nearest_4, a_nearest_4] = compute(p0_nearest_4,d0_nearest_4,p1_nearest_4,d1_nearest_4);

d0_nearest_16 = round((d0)*16)/(16);
d1_nearest_16 = round((d1)*16)/(16);
p0_nearest_16 = round((p0)*16)/(16);
p1_nearest_16 = round((p1)*16)/(16);
[b_nearest_16, a_nearest_16] = compute(p0_nearest_16,d0_nearest_16,p1_nearest_16,d1_nearest_16);

[H_nearest_4, W_nearest_4] = freqz(b_nearest_4,a_nearest_4);
[H_nearest_16, W_nearest_16] = freqz(b_nearest_16,a_nearest_16);

figure;

magH2 = 20*log10(abs(H2));
plot(W2/pi,magH2);

hold on;
magH_nearest_16 = 20*log10(abs(H_nearest_16));
plot(W_nearest_16/pi,magH_nearest_16);

hold on;
magH_nearest_4 = 20*log10(abs(H_nearest_4));
plot(W_nearest_4/pi,magH_nearest_4);

title("Magnitude Response");
xlabel("Normalized Frequency (\times \pi rad/sample)");
ylabel("Magnitude (dB)");
legend("Original", "Nearest 16th", "Nearest 4th");

ax = gca;
ax.YLim = [-50, 10];
ax.XTick = 0:.1:1;

% Part G
figure;
sgtitle("Pole-Zero Plots");

subplot(3,1,1);
zplane(b_computed,a_computed);
title("Original");

subplot(3,1,2);
zplane(b_nearest_16,a_nearest_16);
title("Nearest 16th");

subplot(3,1,3);
zplane(b_nearest_4,a_nearest_4);
title("Nearest 4th");

% The poles are inside the unit for both approximations; we see that the
% parallel allpass realiztion is less sensitive to roundoff error.