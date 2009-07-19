function x = var(s, A)
% VAR  Variance of an observable in a quantum state.
%  x = var(s, A)
%
%  Returns the variance of observable A in the state s.
%  A has to be Hermitian.

% Ville Bergholm 2009


x = ev(s, A^2) - ev(s, A)^2;
