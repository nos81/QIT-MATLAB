function U = rand_SU(n)
% RAND_SU  Random SU(n) matrix.
%  U = rand_SU(n)
%
%  Returns a random special unitary n*n matrix.
%  The matrix is random with respect to the Haar measure.

% Ville Bergholm 2005-2009


U = rand_U(n);
d = det(U)^(1/n); % *exp(i*2*pi*k/n), not unique FIXME
U = U/d;
