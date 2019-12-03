function h= sym2z2d(Hz,zvars,N1N2)
%% function  h= sym2z2d(Hz,zvars,N1N2)
% Convert symbolic 2d z-transform to impulse response matrix
% Hz= 2d z transform
% zvars= symbolic vars in Hz,  [z1,z2]
% N1N2 specifies size:  -N1<=n1<=N1, -N2<=n2<=N2
% h matrix: h[n1,n2]: row index -N1 down to N1, col -N2 to N2
%
% ECE310 DSP  Cooper Union  Fall 2019     Prof FF
z1= zvars(1);
z2= zvars(2);
N1= N1N2(1);
N2= N1N2(2);

M1= 2*N1+1;
M2= 2*N2+1;
h= zeros(M1,M2);
% clear neg powers:
Hz0= simplify(z1^N1*z2^N2*Hz);
Hz1= coeffs(Hz0,z1); %extracts z2-poly as coeff of z1, from z1^0 to z1^max
lenz1= length(Hz1);  % rows lenz1+1 to M1 stay at 0: no work needed
for n1= 1:lenz1
    htmp= sym2poly(Hz1(n1));  % coeffs from z2^max to z2^0, in row vector
    % flip, higher order coeff stay at 0
    h(n1,1:length(htmp))= fliplr(htmp);
end

    
    