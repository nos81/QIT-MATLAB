function [b] = boson_ladder(n)
% BOSON_LADDER  Bosonic ladder operators.
%  b = boson_ladder(n)
%
%  Returns the n-dimensional approximation of the bosonic
%  annihilation operator b for a single bosonic mode in the
%  number basis {|0>, |1>, ..., |n-1>}.
%
%  The corresponding creation operator is b'.

% Ville Bergholm 2009-2010


b = spdiags(sqrt(0:n-1).', 1, n, n);
