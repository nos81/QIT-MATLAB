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

% symmetric
assert_o(fidelity(rho1, rho2), fidelity(rho2, rho1), tol);
assert_o(trace_dist(rho1, rho2), trace_dist(rho2, rho1), tol);

assert_o(fidelity(sigma1, sigma1), 1, tol); % normalized to unity
assert_o(trace_dist(sigma1, sigma1), 0, tol); % distance measure

% unaffected by unitary transformations
assert_o(fidelity(rho1, rho2), fidelity(u_propagate(rho1, U_r), u_propagate(rho2, U_r)), tol);
assert_o(trace_dist(rho1, rho2), trace_dist(u_propagate(rho1, U_r), u_propagate(rho2, U_r)), tol);

% for pure states they're equivalent
assert_o(trace_dist(p1, p2)^2 +fidelity(p1, p2)^2, 1, tol);
% for mixed states, these inequalities hold
assert(sqrt(1-fidelity(rho1, rho2)^2) -trace_dist(rho1, rho2) >= -tol);
assert(1-fidelity(rho1, rho2) -trace_dist(rho1, rho2) <= tol);
% for a pure and a mixed state we get this inequality
assert(1-fidelity(rho1, p1)^2 -trace_dist(rho1, p1) <= tol);


% Test script for entropy.m
% Ville Bergholm 2009

assert_o(entropy(p1), 0, tol); % zero for pure states
assert(entropy(sigma1) >= -tol); % nonnegative

% unaffected by unitary transformations
assert_o(entropy(u_propagate(sigma1, U_s)), entropy(sigma1), tol);


% Test script for ptrace.m
% Ville Bergholm 2008

rho_A = ptrace(rho1, [2]);
% trace of partial trace equals total trace
assert_o(trace(rho1), trace(rho_A), tol)
% partial trace over all subsystems equals total trace
assert_o(trace(rho1), trace(ptrace(rho1, [1 2])), tol)

rho_X = ptrace(sigma1, [2 3 5]);
assert_o(trace(sigma1), trace(rho_X), tol)
assert_o(trace(sigma1), trace(ptrace(sigma1, 1:5)), tol)


% Test script for ptranspose.m
% Ville Bergholm 2008

rho_pt_B = ptranspose(rho1, [2]);
assert_o(trace_dist(rho1, ptranspose(rho_pt_B, [2])), 0, tol)
assert_o(trace(rho1), trace(rho_pt_B), tol)


% Test script for schmidt.m
% Ville Bergholm 2009-2010

[lambda1, u, v] = schmidt(p1, 1);
lambda2 = schmidt(p1, 2);
% squares of schmidt coefficients % sum up to unity
assert_o(norm(lambda1), 1, tol);
% both subdivisions have identical schmidt coefficients
assert_o(norm(lambda1-lambda2), 0, tol);

temp = 0; for k=1:2, temp = temp + kron(lambda1(k)*u(:,k), v(:,k)); end
assert_o(norm(p1.data-temp), 0, tol);

% squared schmidt coefficients equal eigenvalues of partial trace
for k=1:20
  r = normalize(state(rand(30,1)-0.5 +i*(rand(30,1)-0.5), [5 6]));
  x = schmidt(r, [1]).^2;
  temp = ptrace(r, [2]);
  y = sort(eig(temp.data), 'descend');
  assert_o(norm(x-y), 0, tol);
end
