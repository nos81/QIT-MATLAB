% Test script for the harmonic oscillator module.
% Ville Bergholm 2017

tol = qit.tol;

% truncation dimension
m = 100;

% NOTE: Due to the truncation of the Hilbert space, the last
% rows and cols of some operators can be spurious.
drop = 1;
ind = 1:m-drop;
nind = ind(end)+1:m;

I = eye(m);
a = boson_ladder(m);
n = a'*a;
Q = ho.position(m);
P = ho.momentum(m);

% definition of a
temp = Q +1i*P -2*a;
assert_o(norm(temp(ind, ind)), 0, tol);

% ground state
temp = state(0, m);
assert_o(temp.ev(Q^2), 1, tol);
assert_o(temp.ev(P^2), 1, tol);


%% commutation relations

temp = comm(Q,P) -2i*I;
assert_o(norm(temp(ind, ind)), 0, tol);

temp = comm(a,a') -I;
assert_o(norm(temp(ind, ind)), 0, tol);

temp = full(comm(n,a) +a);
assert_o(norm(temp(ind, ind)), 0, tol);


%% position and momentum eigenstates

if 0
q = randn();
p = randn();

psi = ho.position_state(q, m);
psi.ev(Q)-q
psi.ev(Q^2)-q^2

psi = ho.momentum_state(p, m);
psi.ev(P)-p
end


%% coherent states

z = rand_GL(1);
w = rand_GL(1);
D = ho.displace(z, m);
psi = ho.coherent_state(z, m);
xi  = ho.coherent_state(w, m);

% eigenstate of the annihilation op
assert_o(norm(a*psi.data -z*psi.data), 0, tol);

% gaussian overlap
assert_o(xi.data'*psi.data, exp(-1i*imag(w*conj(z))) * exp(-0.5*abs(z-w)^2), tol)

% number op expectation values
mu = psi.ev(n);
assert_o(mu, abs(z)^2, tol);
assert_o(psi.ev(n^2)-mu^2, abs(z)^2, tol);

% coherent states are displaced ground states
temp = state(0, m);
temp = temp.prop(D);
assert_o(norm(psi-temp), 0, tol);
