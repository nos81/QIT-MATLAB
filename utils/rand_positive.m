function A = rand_positive(n)
% RAND_POSITIVE  Random n*n positive semidefinite matrix.
%  A = rand_positive(n)
%
%  Normalized to trace(A) = 1.
%  Since A has all-real eigenvalues, it is Hermitian by construction.

% Ville Bergholm 2008-2012


d = rand_pu(n); % random partition of unity
U = rand_U(n); % random unitary
A = U'*diag(d)*U;

A = (A+A')/2; % eliminate rounding errors
return

% TODO alternative: inverse purification
s = state(0, [n, k]); % rank k state op
s = u_propagate(s, rand_U(n*k)); % expensive and wasteful...
s = ptrace(s, 2);
A = s.data;
end
