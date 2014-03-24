function ret = comm(A, B)
% COMM  Matrix commutator.
%  X = comm(A, B)
%
%  Returns X = AB - BA.

% Ville Bergholm 2014

ret = A*B -B*A;
