% ECE310_Homework_2
% Kevin Kerliu
clear;
close all;
clc;
%%

% Question 1
% Part D

% H
z_H = [1/2; -4];
p_H = [1/3; -5/2; -5/2];
k_H = 1/18;
[b_H,a_H] = zp2tf(z_H,p_H,k_H);
[phi_H,w_H] = phasez(b_H,a_H);
deg_H = phi_H*180/pi;
% A
z_A = [-1/4; -5/2; -5/2];
p_A = [-4; -2/5; -2/5];
k_A = 16/25;
[b_A,a_A] = zp2tf(z_A,p_A,k_A);
[phi_A,w_A] = phasez(b_A,a_A);
deg_A = phi_A*180/pi;
% Hmin
z_Hmin = [1/2; -1/4; 0];
p_Hmin = [1/3; -2/5; -2/5];
k_Hmin = 8/225;
[b_Hmin,a_Hmin] = zp2tf(z_Hmin,p_Hmin,k_Hmin);
[phi_Hmin,w_Hmin] = phasez(b_Hmin,a_Hmin);
deg_Hmin = phi_Hmin*180/pi;

figure;
plot(w_H,deg_H);
hold on;
plot(w_A,deg_A);
hold on;
plot(w_Hmin,deg_Hmin);

title("Phase Responses");
xlabel("Normalized Frequency ( \times \pi rad/sample)");
ylabel("Phase (Degrees)");
legend("H","A","Hmin");