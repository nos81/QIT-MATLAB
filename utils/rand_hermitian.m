function H = rand_hermitian(n)
% RAND_HERMITIAN  Random Hermitian n*n matrix.
%  H = rand_hermitian(n)
%
%  Returns a random Hermitian matrix of size n*n.
%  NOTE: The randomness is not defined in any deeply meaningful sense.

% Ville Bergholm 2008-2011

H = randn(n) +1i*randn(n);
H = H + H'; % make it Hermitian
