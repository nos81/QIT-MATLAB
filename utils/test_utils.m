% Test script for the Quantum Information Toolkit utilities
% Ville Bergholm 2009-2016

tol = qit.tol;

dim = 5;

%% random matrices

H = rand_hermitian(dim);
assert_o(norm(H - H'), 0, tol);

U = rand_U(dim);
assert_o(norm(U*U' -eye(dim)), 0, tol);
assert_o(norm(U'*U -eye(dim)), 0, tol);

U = rand_SU(dim);
assert_o(norm(U*U' -eye(dim)), 0, tol);
assert_o(norm(U'*U -eye(dim)), 0, tol);
assert_o(det(U), 1, tol);

rho = rand_positive(dim);
assert_o(norm(rho-rho'), 0, tol); % hermitian
assert_o(trace(rho), 1, tol); % trace 1
temp = eig(rho);
assert_o(norm(imag(temp)), 0, tol); % real eigenvalues
assert_o(norm(temp-abs(temp)), 0, tol); % nonnegative eigenvalues


%% superoperators

L = rand_U(dim);
R = rand_U(dim);

assert_o(norm(rho -inv_vec(vec(rho))), 0, tol)
assert_o(norm(L*rho*R -inv_vec(lrmul(L, R)*vec(rho))), 0, tol)
assert_o(norm(L*rho -inv_vec(lmul(L)*vec(rho))), 0, tol)
assert_o(norm(rho*R -inv_vec(rmul(R)*vec(rho))), 0, tol)

% vec-superops and Choi matrices should be equivalent
% output and input vector space dims
d_in = 3;
d_out = 4;
L = rand([d_out, d_in].^2);
C = superop_to_choi(L);
temp = gate.copydot(0, 2, d_in);
temp = kron(temp.data, eye(d_out));

rho = rand(d_in);
r1 = inv_vec(L * vec(rho));
r2 = temp' * kron(rho, C) * temp;
assert_o(norm(r1-r2), 0, tol);


%% spectral decomposition

[E, P] = spectral_decomposition(H);
m = length(E); % unique eigenvalues
temp = 0;
for k=1:length(E)
  temp = temp + E(k)*P{k};
end
assert_o(norm(temp-H), 0, tol);


%% tensorsum

A = randn(dim);
B = randn(dim);
C = randn(dim);

assert_o(norm(tensorsum(A, B)'-tensorsum(A', B')), 0, tol);

temp = tensorsum(A, B, C);
assert_o(norm(temp -tensorsum(A, tensorsum(B, C))), 0, tol);
assert_o(norm(temp -tensorsum(tensorsum(A, B), C)), 0, tol);


%% Test script for invariant methods.
% Ville Bergholm 2010

U = rand_U(4); % random two-qubit gate
L = kron(rand_U(2), rand_U(2)); % random local 2-qubit gate
cnot = gate.controlled(qit.sx, 1);
cnot = cnot.data;

% canonical invariants
assert_o(norm(invariant.canonical(L) -[0 0 0]), 0, tol);
assert_o(norm(invariant.canonical(cnot) -[0.5 0 0]), 0, tol);

% Makhlin invariants
c = invariant.canonical(U);
g1 = invariant.makhlin(c);
g2 = invariant.makhlin(U);
assert_o(norm(g1-g2), 0, tol);

% maximum concurrence
assert_o(invariant.max_concurrence(L), 0, tol);
assert_o(invariant.max_concurrence(cnot), 1, tol);
