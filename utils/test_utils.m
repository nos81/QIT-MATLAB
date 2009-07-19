% Test script for the Quantum Information Toolkit utilities
% Ville Bergholm 2009

tol = qit.tol;

dim = 5;


% random matrices
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


% superoperators
L = rand_U(dim);
R = rand_U(dim);

assert_o(norm(rho -inv_vec(vec(rho))), 0, tol)
assert_o(norm(L*rho*R -inv_vec(lrmul(L, R)*vec(rho))), 0, tol)
assert_o(norm(L*rho -inv_vec(lmul(L)*vec(rho))), 0, tol)
assert_o(norm(rho*R -inv_vec(rmul(R)*vec(rho))), 0, tol)

