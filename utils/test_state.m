% Test script for the state class.
% Ville Bergholm 2008-2011

disp('Testing the state class.')

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
dim = [2, 2];
p = state(0, dim);

p1 = prop(p, rand_SU(prod(dim)));
p2 = prop(p, rand_SU(prod(dim)));
U_p = rand_U(prod(dim));

% TODO concurrence, fix_phase, locc_convertible, lognegativity, measure,
% negativity,


%% propagation

dim = [2, 3];
s = state('02', dim);

% with lmaps
q = s.prop(gate.swap(2,3));
assert_o(norm(q-state('20', fliplr(dim))), 0, tol);


%% generalized Bloch vectors

temp = bloch_vector(sigma1, true);
assert_o(norm(bloch_state(temp) -sigma1), 0, tol); % consistency
temp = temp(:);
assert_o(norm(imag(temp)), 0, tol); % correlation tensor is real
assert(sqrt(prod(sigma1.dim)) -norm(temp, 'fro') >= -tol); % purity limit


%% fidelity, trace_dist

% symmetric
assert_o(fidelity(rho1, rho2), fidelity(rho2, rho1), tol);
assert_o(trace_dist(rho1, rho2), trace_dist(rho2, rho1), tol);

assert_o(fidelity(sigma1, sigma1), 1, tol); % normalized to unity
assert_o(trace_dist(sigma1, sigma1), 0, tol); % distance measure

% unaffected by unitary transformations
assert_o(fidelity(rho1, rho2), fidelity(prop(rho1, U_r), prop(rho2, U_r)), tol);
assert_o(trace_dist(rho1, rho2), trace_dist(prop(rho1, U_r), prop(rho2, U_r)), tol);

% for pure states they're equivalent
assert_o(trace_dist(p1, p2)^2 +fidelity(p1, p2)^2, 1, tol);
% for mixed states, these inequalities hold
assert(sqrt(1-fidelity(rho1, rho2)^2) -trace_dist(rho1, rho2) >= -tol);
assert(1-fidelity(rho1, rho2) -trace_dist(rho1, rho2) <= tol);
% for a pure and a mixed state we get this inequality
assert(1-fidelity(rho1, p1)^2 -trace_dist(rho1, p1) <= tol);


%% entropy

assert_o(entropy(p1), 0, tol); % zero for pure states
assert(entropy(sigma1) >= -tol); % nonnegative

% unaffected by unitary transformations
assert_o(entropy(prop(sigma1, U_s)), entropy(sigma1), tol);


%% partial trace

rho_A = ptrace(rho1, [2]);
% trace of partial trace equals total trace
assert_o(trace(rho1), trace(rho_A), tol)
% partial trace over all subsystems equals total trace
assert_o(trace(rho1), trace(ptrace(rho1, [1 2])), tol)

rho_X = ptrace(sigma1, [2 3 5]);
assert_o(trace(sigma1), trace(rho_X), tol)
assert_o(trace(sigma1), trace(ptrace(sigma1, 1:5)), tol)


%% partial transpose

rho_pt_B = ptranspose(rho1, [2]);
assert_o(trace_dist(rho1, ptranspose(rho_pt_B, [2])), 0, tol)
assert_o(trace(rho1), trace(rho_pt_B), tol)


%% schmidt

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


%% reorder

dim = [2 5 1];
A = rand(dim(1));
B = rand(dim(2));
C = rand(dim(3));
T1 = state(mkron(A, B, C), dim);
T2 = reorder(T1, [3 1 2]);
assert_o(norm(mkron(C, A, B) - T2.data), 0, tol);
T2 = reorder(T1, [2 1 3]);
assert_o(norm(mkron(B, A, C) - T2.data), 0, tol);
T2 = reorder(T1, [3 2 1]);
assert_o(norm(mkron(C, B, A) - T2.data), 0, tol);
