function [H_mag, H_phase] = analog_filter(b, a, w, wp_vec, ws_vec, rs, rp, filter)

H = freqs(b, a, w);
H_mag = 20*log10(abs(H));
H_phase = 180/pi*unwrap(angle(H));

figure;
sgtitle(filter);
subplot(2,1,1);
plot(w/(2*pi*1000),H_mag);
title("Magnitude Response");
xlabel("Frequency (kHz)");
ylabel("Magnitude (dB)");
ylim([-100, 2]);
line(wp_vec/(2*pi*1000), [0, 0], 'Color', 'blue', 'LineStyle', '--');
line(wp_vec/(2*pi*1000), [-rp, -rp], 'Color', 'blue', 'LineStyle', '--');
line([0 ws_vec(1)]/(2*pi*1000), [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([ws_vec(2) w(end)]/(2*pi*1000), [-rs,-rs], 'Color', 'blue', 'LineStyle', '--');
line([wp_vec(1) wp_vec(1)]/(2*pi*1000), [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([wp_vec(2) wp_vec(2)]/(2*pi*1000), [-100 2], 'Color', 'blue', 'LineStyle', '--');
line([ws_vec(1) ws_vec(1)]/(2*pi*1000), [-100 2], 'Color', 'blue', 'LineStyle', '--')
line([ws_vec(2) ws_vec(2)]/(2*pi*1000), [-100 2], 'Color', 'blue', 'LineStyle', '--');

subplot(2,1,2);
plot(w/(2*pi*1000),H_phase);
title("Filter Phase Response");
xlabel("Frequency (kHz)");
ylabel("Phase (Degrees)");

end