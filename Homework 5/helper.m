function [gn1n2, Gz1z2] = helper(h_n)
    % h_n is assumed to be a column vector (transpose your input when passing it in if not)
    
    % Part 1
    syms w;
    L = length(h_n); % Length of the filter (Will always be odd; we have defined the support as {-N <= n <= N})
    N = (L-1)/2; % Highest order Chebyshev polynomial
    T = chebyshevT(1:N,cos(w)); % Vector of symbolic Chebyshev polynomials
    
    % Define H(w) 
    H(w) = 2*h_n((L+3)/2:L)'*T'; % Summation of the Chebyshev polynomials and respective coefficients
    H(w) = H(w) + h_n((L+1)/2); % Add the center coefficient (h_0)
    
    % Compute H(w) over [-pi,pi] for 128 equally spaced points
    W = linspace(-pi,pi,128);
    H_freq = H(W);
    H_freq = double(H_freq);

    % Graph the magnitude response of |H(w)| on a linear scale
    figure;
    plot(W,abs(H_freq));
    title("Magnitude Response of |H(w)| on Linear Scale");
    xlabel("Frequency (rad)");
    ylabel("Magnitude (linear)");
    
    
    % Part 2
    syms w1 w2;
    F(w1,w2) = .5*(-1 + cos(w1) + cos(w2) + cos(w1)*cos(w2));
    T = chebyshevT(1:N,F(w1,w2)); % Vector of symbolic Chebyshev polynomials
    
    % Define G(w1,w2) 
    G(w1,w2) = 2*h_n((L+3)/2:L)'*T'; % Summation of the Chebyshev polynomials and respective coefficients
    G(w1,w2) = G(w1,w2) + h_n((L+1)/2); % Add the center coefficient (h_0)
    
    % Compute H(w) over [-pi,pi]^2 for a 128 x 128 grid
    [X,Y] = meshgrid(W,W);
    G_freq = G(X,Y);
    G_freq = double(G_freq);
    
    % Graph the magnitude response of |G(w1,w2)| on a linear scale
    figure;
    plot3(X,Y,abs(G_freq));
    title("Magnitude Response of |G(w1,w2)| on Linear Scale");
    xlabel("Frequency: w1 (rad)");
    ylabel("Frequency: w2 (rad)");
    zlabel("Magnitude (linear)");
    
    
    % Part 3
    syms z1 z2;
    Fz(z1,z2) = .5*(-1 + (z1 + z1^-1)/2 + (z2 + z2^-1)/2 + ((z1 + z1^-1)/2)*((z2 + z2^-1)/2));
    T = chebyshevT(1:N,Fz(z1,z2)); % Vector of symbolic Chebyshev polynomials
    Gz1z2 = h_n((L+1)/2) + 2*sum(h_n((L+3)/2:L).'.*T);
    zvars = [z1 z2];   
    N1N2 = [N,N];
    gn1n2 = sym2z2d(Gz1z2,zvars,N1N2);
%     
end