function f = fidelity(A, B)
% FIDELITY  Relative fidelity of two state operators/matrices.
%  f = fidelity(A, B)
%
%  Returns trace(A*B).

% Ville Bergholm 2007-2009

f = trace(A*B);
