% ECE310_Homework_3
% Kevin Kerliu, Shyam Paidipati
clear;
close all;
clc;
%% Question 1 Part A
% Defining stopband and pasband ranges (rad/sec)
ws1 = 2*pi*1e6;
ws2 = 2*pi*1.6e6;
ws_vec = [ws1 ws2];
wp1 = 2*pi*1.2e6;
wp2 = 2*pi*1.5e6;
wp_vec = [wp1 wp2];

% Defining B(MHz) and w0(rad/sec)
B = (wp2-wp1);
w0 = sqrt(wp1*wp2); 

% Prototype lowpass specs
lpspec_wp = min(abs((wp1^2-w0^2)/(B*wp1)),abs((wp2^2-w0^2)/(B*wp2)));
lpspec_ws = min(abs((ws1^2-w0^2)/(B*ws1)),abs((ws2^2-w0^2)/(B*ws2)));

%% Question 1 Part B
% delta_rs: delta stopband ripple  delta_rp: delta passband ripple
rs = 30;
rp = 2;
delta_rs = 10^(rs/10)-1;
delta_rp = 10^(rp/10)-1;

butterworth_order = ceil(0.5*log10(delta_rs/delta_rp)/log10(lpspec_ws/lpspec_wp));
cheb1_order = ceil(acosh(sqrt(delta_rs/delta_rp))/acosh(lpspec_ws/lpspec_wp));
cheb2_order = cheb1_order;

% Order of prototype lowpass filters
% Butterworth: 9
% Chebyshev I: 5
% Chebyshev II: 5

% Order of bandpass filters
% Butterworth: 18
% Chebyshev I: 10
% Chebyshev II: 10

%% Question 1 Part C
%buttord, butter
[n_butterord,Wn_butterord] = buttord(wp_vec,ws_vec,rp,rs,'s');
[b_butter,a_butter] = butter(n_butterord,Wn_butterord,'bandpass','s');
%cheb1ord, cheby1
[n_cheb1ord,Wn_cheb1ord] = cheb1ord(wp_vec,ws_vec,rp,rs,'s');
[b_cheby1,a_cheby1] = cheby1(n_cheb1ord,rp,wp_vec,'bandpass','s');
%cheb2ord, cheby2
[n_cheb2ord,Wn_cheb2ord] = cheb2ord(wp_vec,ws_vec,rp,rs,'s');
[b_cheby2,a_cheby2] = cheby2(n_cheb2ord,rs,ws_vec,'bandpass','s');
%ellipord, ellip
[n_ellipord,Wn_ellipord] = ellipord(wp_vec,ws_vec,rp,rs,'s');
[b_ellip,a_ellip] = ellip(n_ellipord,rp,rs,wp_vec,'bandpass','s');

%% Question 1 Part D
figure;
zplane(b_butter,a_butter);
title("Bandpass Butterworth");
figure;
zplane(b_cheby1,a_cheby1);
title("Bandpass Chebyshev I");
figure;
zplane(b_cheby2,a_cheby2);
title("Bandpass Chebyshev II");
figure;
zplane(b_ellip,a_ellip);
title("Bandpass Elliptic");

%% Question 1 Part E
[n2_ellipord,Wn2_ellipord] = ellipord(lpspec_wp,lpspec_ws,rp,rs,'s');
[z_lpellip,p_lpellip,~] = ellip(n2_ellipord,rp,rs,lpspec_wp,'low','s');
figure;
zplane(z_lpellip, p_lpellip);
title("Lowpass Elliptic");

%% Question 1 Part F

w = 0:100:2*pi*3e6; 
[H_butter_mag, H_butter_phase] = analog_filter(b_butter, a_butter, w, wp_vec, ws_vec, rs, rp, "Butterworth");
[H_cheby1_mag , H_cheby1_phase] = analog_filter(b_cheby1,a_cheby1, w, wp_vec, ws_vec, rs, rp, "Chebyshev Type I");
[H_cheby2_mag, H_cheby2_phase] = analog_filter(b_cheby2, a_cheby2, w, wp_vec, ws_vec, rs, rp, "Chebyshev Type II");
[H_ellip_mag, H_ellip_phase] = analog_filter(b_ellip, a_ellip, w, wp_vec, ws_vec, rs, rp, "Elliptical");

%% Question 1 Part G

figure;
plot(w/(2*pi*1000),H_butter_mag);
hold on;
plot(w/(2*pi*1000),H_cheby1_mag);
hold on;
plot(w/(2*pi*1000),H_cheby2_mag);
hold on;
plot(w/(2*pi*1000),H_ellip_mag);
line(wp_vec/(2*pi*1000), [0, 0], 'Color', 'blue', 'LineStyle', '--');
line(wp_vec/(2*pi*1000), [-rp, -rp], 'Color', 'blue', 'LineStyle', '--');
xlim(wp_vec/(2*pi*1000));
title("Superimposed Magnitude Responses");
xlabel("Frequency (kHz)");
ylabel("Magnitude (dB)");
legend("Butterworth", "Chebyshev Type I", "Chebyshev Type II", "Elliptical", ...
    "Passband Ripple Upper Limit", "Passband Ripple Lower Limit");

%% Question 1 Part H
ws1_index = find(w < ws_vec(1));
ws2_index = find(w > ws_vec(2), 1);
butter_atten_ws1 = max(H_butter_mag(ws1_index));
butter_atten_ws2 = max(H_butter_mag(ws2_index));
cheby1_atten_ws1 = max(H_cheby1_mag(ws1_index));
cheby1_atten_ws2 = max(H_cheby1_mag(ws2_index));
cheby2_atten_ws1 = max(H_cheby2_mag(ws1_index));
cheby2_atten_ws2 = max(H_cheby2_mag(ws2_index));
ellip_atten_ws1 = max(H_ellip_mag(ws1_index));
ellip_atten_ws2 = max(H_ellip_mag(ws2_index));
table(butter_atten_ws1, cheby1_atten_ws1, cheby2_atten_ws1, ellip_atten_ws1)
table(butter_atten_ws2, cheby1_atten_ws2, cheby2_atten_ws2, ellip_atten_ws2)

%% Question 2 Part A

% Sampling frequency, stopband and passband ripple
fs = 6e6;
rs = 30;
rp = 2;

% Converting specifications 
ws1_d = (1e6)/(fs/2);
ws2_d = (1.6e6)/(fs/2);
ws_vec_d = [ws1_d ws2_d];
wp1_d = (1.2e6)/(fs/2);
wp2_d = (1.5e6)/(fs/2);
wp_vec_d = [wp1_d wp2_d];

% Defining B(MHz) and w0(rad/sec)
B_d = (wp2_d-wp1_d);
w0_d = sqrt(wp1_d*wp2_d); 

% Defining lowpass specs
lpspec_wp_d = min(abs((wp1_d^2-w0_d^2)/(B_d*wp1_d)),abs((wp2_d^2-w0_d^2)/(B_d*wp2_d)));
lpspec_ws1_d = abs((ws1_d^2-w0_d^2)/(B_d*ws1_d));
lpspec_ws2_d = abs((ws2_d^2-w0_d^2)/(B_d*ws2_d));
lpspec_ws_d = min(lpspec_ws1_d,lpspec_ws2_d);

%% Question 2 Part B

%buttord, butter
[n_butterord_d,Wn_butterord_d] = buttord(wp_vec_d,ws_vec_d,rp,rs);
[z_butter_d,p_butter_d,k_butter_d] = butter(n_butterord_d,Wn_butterord_d,'bandpass');
%cheb1ord, cheby1
[n_cheb1ord_d,Wn_cheb1ord_d] = cheb1ord(wp_vec_d,ws_vec_d,rp,rs);
[z_cheby1_d,p_cheby1_d,k_cheby1_d] = cheby1(n_cheb1ord_d,rp,wp_vec_d,'bandpass');
%cheb2ord, cheby2
[n_cheb2ord_d,Wn_cheb2ord_d] = cheb2ord(wp_vec_d,ws_vec_d,rp,rs);
[z_cheby2_d,p_cheby2_d,k_cheby2_d] = cheby2(n_cheb2ord_d,rs,ws_vec_d,'bandpass');
%ellipord, ellip
[n_ellipord_d,Wn_ellipord_d] = ellipord(wp_vec_d,ws_vec_d,rp,rs);
[z_ellip_d,p_ellip_d,k_ellip_d] = ellip(n_ellipord_d,rp,rs,wp_vec_d,'bandpass');

% Order of the prototype filters
% Butterworth: 8
% Chebyshev I: 5
% Chebyshev II: 5
% Elliptic: 3

% Order of the final filters
% Butterworth: 16
% Chebyshev I: 10
% Chebyshev II: 10
% Elliptic: 6

%% Question 2 Part C

fs1 = 1;
fs2 = 1.6;
fs_vec = [fs1 fs2];
fp1 = 1.2;
fp2 = 1.5;
fp_vec = [fp1 fp2];
n = 1e6;

[H_butter_d, f_butter_d] = digital_filter(z_butter_d, p_butter_d, k_butter_d, n, fs_vec, fp_vec, rs, rp, fs, "Butterworth");

[H_cheby1_d, f_cheby1_d] = digital_filter(z_cheby1_d, p_cheby1_d, k_cheby1_d, n, fs_vec, fp_vec, rs, rp, fs, "Chebyshev Type I");

[H_cheby2_d, f_cheby2_d] = digital_filter(z_cheby2_d, p_cheby2_d, k_cheby2_d, n, fs_vec, fp_vec, rs, rp, fs, "Chebyshev Type II");

[H_ellip_d, f_ellip_d] = digital_filter(z_ellip_d, p_ellip_d, k_ellip_d, n, fs_vec, fp_vec, rs, rp, fs, "Elliptic");

%% Question 2 Part D

figure;
sgtitle("Pole-Zero Plots");

subplot(2, 2, 1);
zplane(z_butter_d, p_butter_d);
title("Butterworth");

subplot(2, 2, 2);
zplane(z_cheby1_d, p_cheby1_d);
title("Chebyshev Type I");

subplot(2, 2, 3);
zplane(z_cheby2_d, p_cheby2_d);
title("Chebyshev Type II");

subplot(2, 2, 4);
zplane(z_ellip_d, p_ellip_d);
title("Elliptic");

%% Question 3 Part A

% Using b_cheby1 and a_cheby from Question 1 and fs from Question 2
[bz_cheby1,az_cheby1] = impinvar(b_cheby1,a_cheby1, fs);
[bz_cheby2,az_cheby2] = impinvar(b_cheby2,a_cheby2, fs);
sizeb_1 = size(bz_cheby1)
sizea_1 = size(az_cheby1)
sizeb_2 = size(bz_cheby2)
sizea_2 = size(az_cheby2)

% sizeb_1 is 11 
% sizea_1 is 11
% sizeb_2 is 11
% sizea_2 is 11

% The lengths of the vectors are all the same; they all have a length of 11
% This corresponds to a filter order of 10, which is the same order as the
% analog filters.

%% Question 3 Part B

n = 1e6;
[Hz_cheby1,fz_cheby1] = freqz(bz_cheby1, az_cheby1,n,fs/1e6);
[Hz_cheby2,fz_cheby2] = freqz(bz_cheby2, az_cheby2,n,fs/1e6);
Hz_cheby1_mag = 20*log10(abs(Hz_cheby1));
Hz_cheby2_mag = 20*log10(abs(Hz_cheby2));
Hz_cheby1_phase = 180/pi*unwrap(angle(Hz_cheby1));
Hz_cheby2_phase = 180/pi*unwrap(angle(Hz_cheby2));

% Magnitude Responses
figure;
sgtitle("Magnitude Responses");
subplot(2,1,1);
plot(fz_cheby1,Hz_cheby1_mag);
title("Chebyshev Type I");
xlabel("Frequency (MHz)");
ylabel("Magnitude (dB)");
ylim([-100, 2]);
line(fp_vec, [0, 0], 'Color', 'blue', 'LineStyle', '--');
line(fp_vec, [-rp, -rp], 'Color', 'blue', 'LineStyle', '--');
line([0 fs_vec(1)], [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([fs_vec(2) fz_cheby1(end)], [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([fp_vec(1) fp_vec(1)], [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([fp_vec(2) fp_vec(2)], [-100 2], 'Color', 'blue', 'LineStyle', '--');
line([fs_vec(1) fs_vec(1)], [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([fs_vec(2) fs_vec(2)], [-100 2], 'Color', 'blue', 'LineStyle', '--');

subplot(2,1,2);
plot(fz_cheby2,Hz_cheby2_mag);
title("Chebyshev Type II");
xlabel("Frequency (MHz)");
ylabel("Magnitude (dB)");
ylim([-50, 2]);
line(fp_vec, [0, 0], 'Color', 'blue', 'LineStyle', '--');
line(fp_vec, [-rp, -rp], 'Color', 'blue', 'LineStyle', '--');
line([0 fs_vec(1)], [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([fs_vec(2) fz_cheby2(end)], [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([fp_vec(1) fp_vec(1)], [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([fp_vec(2) fp_vec(2)], [-100 2], 'Color', 'blue', 'LineStyle', '--');
line([fs_vec(1) fs_vec(1)], [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([fs_vec(2) fs_vec(2)], [-100 2], 'Color', 'blue', 'LineStyle', '--');

% Both magnitude responses are NOT equiripple in the passband

% Phase Responses
figure;
sgtitle("Phase Responses");
subplot(2,1,1);
plot(fz_cheby1,Hz_cheby1_phase);
title("Chebyshev Type I");
xlabel("Frequency (MHz)");
ylabel("Phase (Degrees)");
subplot(2,1,2);
plot(fz_cheby2,Hz_cheby2_phase);
title("Chebyshev Type II");
xlabel("Frequency (MHz)");
ylabel("Phase (Degrees)");

%Checking attenuation 
% Cheby Type I
% Stopband
index_s1_c1 = find(fz_cheby1 >= fs1, 1);
rs_s1_c1 = abs(max(Hz_cheby1_mag(1:index_s1_c1)))
index_s2_c1 = find(fz_cheby1 >= fs2, 1);
rs_s2_c1 = abs(max(Hz_cheby1_mag(index_s2_c1:end)))
% Passband
index_p1_c1 = find(fz_cheby1 >= fp1, 1);
rp_p1_c1 = abs(Hz_cheby1_mag(index_p1_c1))
index_p2_c1 = find(fz_cheby1 >= fp2, 1);
rp_p2_c1 = abs(Hz_cheby1_mag(index_p2_c1))

% Cheby Type 2
% Stopband
index_s1_c2 = find(fz_cheby2 >= fs1, 1);
rs_s1_c2 = abs(max(Hz_cheby2_mag(1:index_s1_c2)))
index_s2_c2 = find(fz_cheby2 >= fs2, 1);
rs_s2_c2 = abs(max(Hz_cheby2_mag(index_s2_c2:end)))
% Passband
index_p1_c2 = find(fz_cheby2 >= fp1, 1);
rp_p1_c2 = abs(Hz_cheby2_mag(index_p1_c2))
index_p2_c2 = find(fz_cheby2 >= fp2, 1);
rp_p2_c2 = abs(Hz_cheby2_mag(index_p2_c2))

% The Chebyshev Type II meets none of the specifications!

%% Question 3 Part C

figure;
sgtitle("Pole-Zero Plots");
subplot(2,1,1);
zplane(bz_cheby1,az_cheby1);
title("Chebyshev Type I");
subplot(2,1,2);
zplane(bz_cheby2,bz_cheby2);
title("Chebyshev Type II");

% Chebyshev Type I
% The poles remain inside the unit circle
% The zeros moved, some are near original spots

% Chebyshev Type II
% Some of the poles are slgihtly outside the unit circle
% The poles and zeros ovelap!
% The zeros moved, some are near original spots, one is very far away

%% Question 3 Part D

[A_c1,B_c1,C_c1,D_c1] = tf2ss(b_cheby1,a_cheby1);
[A_c2,B_c2,C_c2,D_c2] = tf2ss(b_cheby2,a_cheby2);
sys_c1 = ss(A_c1,B_c1,C_c1,D_c1);
sys_c2 = ss(A_c2,B_c2,C_c2,D_c2);

figure;
sgtitle("Chebyshev Type I")
subplot(2,1,1);
impulse(sys_c1);
title("Analog Impulse Response");
subplot(2,1,2);
impz(bz_cheby1,az_cheby1,450,6e6);
title("Digital Impulse Response");

figure;
sgtitle("Chebyshev Type II");
subplot(2,1,1);
impulse(sys_c2);
title("Analog Impulse Response");
subplot(2,1,2);
impz(bz_cheby2,az_cheby2,180,6e6);
title("Digital Impulse Response");

%% Question 4 Part A

L = 31;
r = 30;
cheby_w = chebwin(L,r);
cheby_w_normalized = cheby_w/sum(cheby_w);

%% Question 4 Part B

cheby_z = roots(cheby_w_normalized);
figure;
zplane(cheby_z);
title("Chebyshev Window Zeros");

%% Question 4 Part C

angle_cheby_w = angle(cheby_z);
width = 2*min(abs(angle_cheby_w));

%% Question 4 Part D

N = 1000;
cheby_w_fft = fftshift(fft(cheby_w_normalized, N));
cheby_w_mag = 20*log10(abs(cheby_w_fft));
figure;
plot(cheby_w_mag);
xlabel("Radians (rad)");
ylabel("Magnitude (dB)");
title("Chebyshev Window Magnitude Response");

%% Question 4 Part E

beta = 3.15;
kaiser_w = kaiser(L, beta);
kaiser_w_normalized = kaiser_w/sum(kaiser_w);

kaiser_z = roots(kaiser_w_normalized);
figure;
zplane(kaiser_z);
title("Kaiser Window Zeros");

angle_kaiser_w = angle(kaiser_z);
width = 2*min(abs(angle_kaiser_w));

kaiser_w_fft = fftshift(fft(kaiser_w_normalized, N));
kaiser_w_mag = 20*log10(abs(kaiser_w_fft));
figure;
plot(kaiser_w_mag);
xlabel("Radians (rad)");
ylabel("Magnitude (dB)"); 
title("Kaiser Window Magnitude Response");

%% Question 4 Part F

figure;
plot(cheby_w_mag);
hold on;
plot(kaiser_w_mag);
xlabel("Radians (rad)");
ylabel("Magnitude (dB)");
title("Chebyshev Window Magnitude Response");

figure;
zplane(cheby_z,kaiser_z);
title("Chebyshev and Kaiser Window Zeros");
legend("Chebyshev","Kaiser");


%% Question 4 Part G

beta; % 3.15

%% Question 4 Part H

[kaiser_peaks, kaiser_loc] = findpeaks(kaiser_w_mag);
[kaiser_peak, kaiser_index] = max(kaiser_peaks);
kaiser_peak_sidelobe_level = kaiser_peaks(kaiser_index + 1)

[cheby_peaks, cheby_loc] = findpeaks(cheby_w_mag);
[cheby_peak, cheby_index] = max(cheby_peaks);
cheby_peak_sidelobe_level = cheby_peaks(cheby_index + 1)

%% Question 4 Part I

% Chebyshev Window
abs_cheby_fft = abs(cheby_w_fft);
cheby_index_max_peak = find(abs_cheby_fft(cheby_loc) == cheby_peak);
cheby_loc_min = islocalmin(abs_cheby_fft);
cheby_right_null = find(cheby_loc_min(cheby_loc(cheby_index):end) == 1 , 1) + cheby_loc(cheby_index);
cheb_left_null = find(cheby_loc_min(1:cheby_loc(cheby_index)) == 1);
cheb_left_null = cheb_left_null(end);
% Computing energies
cheby_tot_energy = sum(abs_cheby_fft.^2);
cheby_main_lobe_energy = sum(abs_cheby_fft(cheb_left_null:cheby_right_null).^2);
cheby_side_lobe_fraction_energy = (cheby_tot_energy - cheby_main_lobe_energy) / cheby_tot_energy

% Kaiser Window
abs_kaiser_fft = abs(kaiser_w_fft);
kaiser_index_max_peak = find(abs_kaiser_fft(kaiser_loc) == kaiser_peak);
kaiser_loc_min = islocalmin(abs_kaiser_fft);
kaiser_right_null = find(kaiser_loc_min(kaiser_loc(kaiser_index):end) == 1 , 1) + kaiser_loc(kaiser_index);
kaiser_left_null = find(kaiser_loc_min(1:kaiser_loc(kaiser_index)) == 1);
kaiser_left_null = kaiser_left_null(end);
% Computing energies
kaiser_tot_energy = sum(abs_kaiser_fft.^2);
kaiser_main_lobe_energy = sum(abs_kaiser_fft(kaiser_left_null:kaiser_right_null).^2);
kaiser_side_lobe_fraction_energy = (kaiser_tot_energy - kaiser_main_lobe_energy) / kaiser_tot_energy

%% Question 5 Part A

fs = 6e6;
rs = 30;
rp = 2;

delta_p = (10^(rp/20)-1) / (10^(rp/20)+1);
delta_s = 10^(-rs/20);

mags = [0 1 0];
f_cutoff = [1e6 1.2e6 1.5e6 1.6e6];
deviation = [delta_s delta_p delta_s];

% Kaiser
[n_k, Wn_k, beta_k, ftype_k] = kaiserord(f_cutoff, mags, deviation, fs);
N = n_k; % save for later
n_k = n_k + rem(n_k,2);
kaiser_b = fir1(n_k, Wn_k, ftype_k, kaiser(n_k+1,beta_k), 'noscale');

%Parks-McClellan
[n_pm,fo_pm,ao_pm,w_pm] = firpmord(f_cutoff,mags,deviation,fs);
parks_b = firpm(n_pm,fo_pm,ao_pm,w_pm);

%% Question 5 Part B

[kaiser_ht,kaiser_t] = impz(kaiser_b);
[parks_ht,parks_t] = impz(parks_b);

figure;
stem(kaiser_ht);
title("Kaiser");

figure;
stem(parks_ht);
title("Parks-McClellan");

[kaiser_Hf,kaiser_f] = freqz(kaiser_b,1,1024,fs);
[parks_Hf,parks_f] = freqz(parks_b,1,1024,fs);

kaiser_H_mag = 20*log10(abs(kaiser_Hf));
parks_H_mag = 20*log10(abs(parks_Hf));

stopband1 = linspace(0,1);
passband = linspace(1.2,1.5);
stopband2 = linspace(1.6,3);

figure;
sgtitle("Magnitude Responses");

subplot(2,1,1);
plot(kaiser_f/1e6,kaiser_H_mag);
hold on;
title("Kaiser");
xlabel("Frequency (MHz)");
ylabel("Magnitude (dB)");
ylim([-50 2]);
plot(stopband1, -rs*ones(length(stopband1)), "k--");
plot(stopband2, -rs*ones(length(stopband2)), "k--");
plot(passband, (max(kaiser_H_mag))*ones(size(passband)), "k--");
plot(passband,(-rp)*(max(parks_H_mag))*ones(size(passband)), "k--");

subplot(2,1,2);
plot(parks_f/1e6,parks_H_mag);
hold on;
plot(stopband1, -rs*ones(length(stopband1)), "k--");
plot(stopband2, -rs*ones(length(stopband2)), "k--");
plot(passband,(max(parks_H_mag))*ones(size(passband)), "k--");
plot(passband,(-rp)*(max(parks_H_mag))*ones(size(passband)), "k--");
ylim([-50 2]);
title("Parks-McClellan");
xlabel("Frequency (MHz)");
ylabel("Magnitude (dB)");

%% Question 5 Part C

kaiser_passband = kaiser_Hf(kaiser_f >= 1.2e6 & kaiser_f <= 1.5e6);
kaiser_f_min = min(kaiser_passband);
kaiser_f_max = max(kaiser_passband);
kaiser_variation = 20*log10(abs(kaiser_f_max)/abs(kaiser_f_min))
kaiser_stopband = kaiser_H_mag(kaiser_f < 1e6 | kaiser_f > 1.6e6);
kaiser_attenuation = abs(max(kaiser_stopband))

parks_passband = parks_Hf(parks_f >= 1.2e6 & parks_f <= 1.5e6);
parks_f_min = min(parks_passband);
parks_f_max = max(parks_passband);
parks_variation = 20*log10(abs(parks_f_max)/abs(parks_f_min))
parks_stopband = parks_H_mag(parks_f < 1e6 | parks_f > 1.6e6);
parks_attenuation = abs(max(parks_stopband))

% Kaiser
% Variation: 0.4785 
% Attenuation: 28.9072

% Parks-McClellan
% Variation: 2.1889, not within spec.
% Attenuation: 29.2078

%% Question 5 Part D

weight = max(delta_p/delta_s, delta_s/delta_p);
checkweight = (weight == w(find(w, 1, "first")));

% Check is good!

%% Question 5 Part E

new_kaiser_b = fir1(N+2,Wn_k,ftype_k,kaiser(3+N,beta_k),"noscale");
[new_kaiser_H,new_kaiser_f] = freqz(new_kaiser_b,1,512,fs);
new_kaiser_H_mag = 20*log10(abs(new_kaiser_H));

[new_n_pm,new_fo_pm,new_ao_pm,new_w_pm] = firpmord(f_cutoff, mags, deviation, fs);
new_parks_b = firpm(new_n_pm + 2,new_fo_pm,new_ao_pm,new_w_pm);
[new_parks_H,new_parks_f] = freqz(new_parks_b,1,1024,fs);
new_parks_H_mag = 20*log10(abs(new_parks_H));

figure;
plot(new_kaiser_f/1e6,new_kaiser_H_mag);
title("Kaiser Magnitude Response");
ylim([-100 2]);
xlabel("Frequency (MHz)");
ylabel("Magnitude (dB)");
line(fp_vec, [0, 0], 'Color', 'blue', 'LineStyle', '--');
line(fp_vec, [-rp, -rp], 'Color', 'blue', 'LineStyle', '--');
line([0 fs_vec(1)], [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([fs_vec(2) new_kaiser_f(end)/1e6], [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([fp_vec(1) fp_vec(1)], [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([fp_vec(2) fp_vec(2)], [-100 2], 'Color', 'blue', 'LineStyle', '--');
line([fs_vec(1) fs_vec(1)], [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([fs_vec(2) fs_vec(2)], [-100 2], 'Color', 'blue', 'LineStyle', '--');

figure;
plot(new_parks_f/1e6,new_parks_H_mag);
title("Park-McClellan Magnitude Response");
ylim([-100 2]);
xlabel("Frequency (MHz)");
ylabel("Magnitude (dB)");
line(fp_vec, [rp, rp], 'Color', 'blue', 'LineStyle', '--');
line(fp_vec, [0, 0], 'Color', 'blue', 'LineStyle', '--');
line(fp_vec, [-rp, -rp], 'Color', 'blue', 'LineStyle', '--');
line([0 fs_vec(1)], [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([fs_vec(2) new_parks_f(end)/1e6], [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([fp_vec(1) fp_vec(1)], [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([fp_vec(2) fp_vec(2)], [-100 2], 'Color', 'blue', 'LineStyle', '--');
line([fs_vec(1) fs_vec(1)], [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([fs_vec(2) fs_vec(2)], [-100 2], 'Color', 'blue', 'LineStyle', '--');

