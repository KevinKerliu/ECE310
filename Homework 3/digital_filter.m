function[H, f] = digital_filter(z, p, k, n, fs_vec, fp_vec, rs, rp, fs, filter);

[b, a] = zp2tf(z, p, k);
[H, f] = freqz(b, a, n, fs);
H_mag = 20*log10(abs(H));
H_phase = 180/pi*unwrap(angle(H));

figure;
sgtitle(filter);
subplot(2,1,1);
plot(f/1e6,H_mag);
title("Magnitude Response");
xlabel("Frequency (MHz)");
ylabel("Magnitude (dB)");
ylim([-100, 2]);
line(fp_vec, [0, 0], 'Color', 'blue', 'LineStyle', '--');
line(fp_vec, [-rp, -rp], 'Color', 'blue', 'LineStyle', '--');
line([0 fs_vec(1)], [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([fs_vec(2) f(end)/1e6], [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([fp_vec(1) fp_vec(1)], [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([fp_vec(2) fp_vec(2)], [-100 2], 'Color', 'blue', 'LineStyle', '--');
line([fs_vec(1) fs_vec(1)], [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([fs_vec(2) fs_vec(2)], [-100 2], 'Color', 'blue', 'LineStyle', '--');

subplot(2,1,2);
plot(f/1e6,H_phase);
title("Filter Phase Response");
xlabel("Frequency (MHz)");
ylabel("Phase (Degrees)");

end