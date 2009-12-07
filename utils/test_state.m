% Test script for the state class.
% Ville Bergholm 2008-2009

tol = qit.tol;

% mixed states
dim = [2 2];
rho1 = state(rand_positive(prod(dim)), dim);
rho2 = state(rand_positive(prod(dim)), dim);
U_r = rand_U(prod(dim));

dim = [2 3 5 2];
sigma1 = state(rand_positive(prod(dim)), dim);
U_s = rand_U(prod(dim));

% pure states
dim = [2 2];
p = state(0, dim);

p1 = u_propagate(p, rand_SU(prod(dim)));
p2 = u_propagate(p, rand_SU(prod(dim)));
U_p = rand_U(prod(dim));



% Test scripts for fidelity.m, trace_dist.m
% Ville Bergholm 2009
assert_o(fidelity(rho1, rho2), fidelity(rho2, rho1), tol);
assert_o(fidelity(sigma1, sigma1), 1, tol);
assert_o(trace_dist(p1, p2)^2 +fidelity(p1, p2)^2, 1, tol);
assert(1-fidelity(rho1, rho2) -trace_dist(rho1, rho2) <= tol);
assert(sqrt(1-fidelity(rho1, rho2)^2) -trace_dist(rho1, rho2) >= -tol);
assert(1-fidelity(rho1, p1)^2 -trace_dist(rho1, p1) <= tol);


% Test script for entropy.m
% Ville Bergholm 2009

assert_o(entropy(p1), 0, tol);
assert(entropy(sigma1) >= -tol);
assert_o(entropy(u_propagate(sigma1, U_s)), entropy(sigma1), tol);


% Test script for ptrace.m
% Ville Bergholm 2008

rho_A = ptrace(rho1, [2]);
assert_o(trace(rho1), trace(rho_A), tol)
assert_o(trace(rho1), trace(ptrace(rho1, [1 2])), tol)

rho_X = ptrace(sigma1, [2 3 5]);
assert_o(trace(sigma1), trace(rho_X), tol)
assert_o(trace(sigma1), trace(ptrace(sigma1, 1:5)), tol)


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
