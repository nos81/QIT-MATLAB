% Test script for Born-Markov methods.

% Ville Bergholm 2009-2010

tol = qit.tol;

dim = 6;

H = rand_hermitian(dim);
D = {rand_hermitian(dim)/10, rand_hermitian(dim)/10};
baths = {markov.bath('ohmic', 1e9, 0.02), markov.bath('ohmic', 1e9, 0.03)};


% jump operators
[dH, A] = markov.ops(H, D);

assert_o(dH(1), 0, tol);
for n=1:length(D)
  temp = A{n,1}; % dH(1) == 0
  for k=2:length(dH)
    temp = temp +A{n,k} +A{n,k}'; % A(-omega) == A'(omega)
  end
  assert_o(norm(temp - D{n}), 0, tol); % Lindblad ops should sum to D
end


% equivalence of Lindblad operators and the Liouvillian superoperator
[LL, H_LS] = markov.lindblad_ops(H, D, baths);
S1 = superop_lindblad(LL, H + H_LS);
S2 = markov.superop(H, D, baths);
assert_o(norm(S1-S2), 0, tol);
