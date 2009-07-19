% Test script for the state class.
% Ville Bergholm 2008-2009

tol = qit.tol;

% mixed states
dim = [2 2];
rho1 = state(rand_positive(prod(dim)), dim);
U1 = rand_U(prod(dim));

dim = [2 3 5 2];
rho2 = state(rand_positive(prod(dim)), dim);
U2 = rand_U(prod(dim));

% pure states
dim = [2 3];
p1 = state(0, dim);
p1.data = rand_SU(prod(dim))*p1.data;



% Test script for entropy.m
% Ville Bergholm 2009

assert_o(entropy(p1), 0, tol);
assert(entropy(rho2) >= 0);
assert_o(entropy(u_propagate(rho2, U2)), entropy(rho2), tol);


% Test script for ptrace.m
% Ville Bergholm 2008

rho_A = ptrace(rho1, [2]);
assert_o(trace(rho1), trace(rho_A), tol)
assert_o(trace(rho1), trace(ptrace(rho1, [1 2])), tol)

rho_X = ptrace(rho2, [2 3 5]);
assert_o(trace(rho2), trace(rho_X), tol)
assert_o(trace(rho2), trace(ptrace(rho2, 1:5)), tol)


% Test script for ptranspose.m
% Ville Bergholm 2008

rho_pt_B = ptranspose(rho1, [2]);
%assert(norm(rho1 - ptranspose(rho_pt_B, [2])) <= tol)
assert_o(trace(rho1), trace(rho_pt_B), tol)


% Test script for schmidt.m
% Ville Bergholm 2009

lambda1 = schmidt(p1, 1);
lambda2 = schmidt(p1, 2);
assert_o(norm(lambda1), 1, tol)
assert_o(norm(lambda1-lambda2), 0, tol);
