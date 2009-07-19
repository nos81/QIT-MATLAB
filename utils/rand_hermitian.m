function H = rand_hermitian(n)
% RAND_HERMITIAN  Random Hermitian n*n matrix.
%  H = rand_hermitian(n)
%
%  Returns a random Hermitian matrix of size n*n.
%  NOTE: The randomness is not defined in any deeply meaningful sense.

% Ville Bergholm 2008-2009

H = (rand(n) - 0.5) +i*(rand(n) - 0.5);
H = H + H'; % make it Hermitian
