% ECE310_Homework_5
% Jungang Leon Fang, Kevin Kerliu, Shyam Paidipati

close all;
clear;
clc;

%% Part A

% 1
F = @(w1,w2) (1/2)*(-1+cos(w1)+cos(w2)+ cos(w1).*cos(w2));
W = linspace(-pi,pi,128);
test1 = F(W,W);
test2 = F(-W,-W);
test3 = F(W,-W);
test4 = F(-W,W);
% Check if test1 == test 2 == test3 == test4.
% They are in fact the same; so we have verified the symmetry.

% 2
test5 = F(pi,W);
test6 = F(W,pi);
% w = pi
% cos(w) = -1
% Check if test5 == test6 == cos(pi).
% They are in fact the same; so we have verified that either 
% w1 = pi or w2 = pi corresponds to w = pi.

%3 
% Show F(w1,0) = cos(w1) and F(0,w2) = cos(w2)
% F(w1,0) = 0.5*(-1+cos(w1)+cos(0)+cos(w1)*cos(0))
%         = 0.5*(-1+cos(w1)+1+cos(w1)*1))
%         = cos(w1)
% F(0,w2) = 0.5*(-1+cos(0)+cos(w2)+cos(0)*cos(w2))
%         = 0.5*(-1+1+cos(w2)+1*cos(w2))
%         = cos(w2)

% Show G(w1,0) = H(w1) and G(0,w2) = H(w2)
% G(w1,0) = h0 + 2*sum(m=1:N)(h_m*T_m(F(w1,0))
%         = h0 + 2*sum(m=1:N)(h_m*T_m(cos(w1))
%         = H(w1)
% G(0,w2) = h0 + 2*sum(m=1:N)(h_m*T_m(F(0,w2))
%         = h0 + 2*sum(m=1:N)(h_m*T_m(cos(w2))
%         = H(w2)
%% Part B

v = linspace(-pi,pi,128);
[X,Y] = meshgrid(v,v);
M = F(X,Y);
M = double(M);
figure;
contour(X,Y,M);
title("Contour plot of F(w1,w2)");
% checked that all values are bouned between -1 and 1

%% Part C

syms z1 z2
F_z(z1,z2) = (1/2)*(-1+0.5*(z1 + z1^(-1))+0.5*(z2 + z2^(-1))+ (0.5*(z2 + z2^(-1)))*(0.5*(z1 + z1^(-1))));

%% Part D

% The support is {-N <= n1 >= N} for the first dimension and 
%                {-N <= n2 <= N} for the second dimension.
% We can also say the support is (2N + 1) x (2N + 1). 

%% Part E

% see helper function

%% Part F

N = 7;
hamming_window = hamming(N);
[gn1n2,Gz1z2] = helper(hamming_window);

%% Part G
% Construct brick wall (ideal) filter
index = (-3:1:3);
h_0 = (sin(pi*index*2*(1/3)))./(index*pi) - (sin((pi/3)*index))./(pi*index);
h_0(4) = 1/3;
trunc = (h_0').*(hamming(7));
[g,G] = helper(trunc);