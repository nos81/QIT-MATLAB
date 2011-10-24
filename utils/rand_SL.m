function S = rand_SL(n)
% RAND_SL  Random SL(n,C) matrix.
%  S = rand_SL(n)
%
%  Returns a random special linear n*n matrix.
%  NOTE: The randomness is not defined in any deeply meaningful sense.

% Ville Bergholm 2011


S = randn(n) +1i*randn(n);
S = S / det(S)^(1/n);
