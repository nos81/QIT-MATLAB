function ret = acomm(A, B)
% ACOMM  Matrix anticommutator.
%  X = acomm(A, B)
%
%  Returns X = AB + BA.

% Ville Bergholm 2014

ret = A*B +B*A;
