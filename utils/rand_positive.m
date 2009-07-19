function A = rand_positive(n)
% RAND_POSITIVE  Random n*n positive semidefinite matrix.
%  A = rand_positive(n)
%
%  Normalized as trace(A) = 1.
%  Since the matrix has purely real eigenvalues, it is also
%  Hermitian by construction.

% Ville Bergholm 2008-2009


p = sort(rand(n-1,1));  % n-1 points in [0,1]
d = sort([p;1]-[0;p]);  % n deltas between points = partition of unity

U = rand_U(n); % random unitary
A = U'*diag(d)*U;

A = (A+A')/2; % eliminate rounding errors
