% Test script for Born-Markov methods.

% Ville Bergholm 2009-2017

tol = qit.tol;

dim = 6;

TU = 1e-9; % s

H = rand_hermitian(dim);
D = {rand_hermitian(dim)/10, rand_hermitian(dim)/10};
baths = {markov.bath('ohmic', 'boson', TU, 0.02),...
         markov.bath('ohmic', 'fermion', TU, 0.03)};


% equivalence of Lindblad operators and the Liouvillian superoperator
[LL, H_LS] = markov.lindblad_ops(H, D, baths);
S1 = superop_lindblad(LL, H + H_LS);
S2 = markov.superop(H, D, baths);
assert_o(norm(S1-S2), 0, tol);
