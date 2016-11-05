function S = rand_GL(n)
% RAND_GL  Random GL(n,C) matrix.
%  S = rand_GL(n)
%
%  Returns a random complex general linear n*n matrix.
%  NOTE: The randomness is not defined in any deeply meaningful sense.
%  The det(S) ~= 0 condition is for now fulfilled only statistically,
%  i.e. you almost never obtain a non-invertible matrix.

% Ville Bergholm 2016


S = randn(n) +1i*randn(n);
