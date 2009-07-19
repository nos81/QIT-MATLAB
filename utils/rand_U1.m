function U = rand_U1(n)
% RAND_U1  Random diagonal unitary matrix.
%  U = rand_U1(n)
%
%  Returns a random diagonal unitary n*n matrix.
%  The matrix is random with respect to the Haar measure.

% Ville Bergholm 2005-2009


U = diag(exp(i*2*pi*rand(n,1)));
