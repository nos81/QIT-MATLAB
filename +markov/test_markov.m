% Test script for Markov methods.

% Ville Bergholm 2009-2010

tol = qit.tol;

dim = 5;


% Lindblad operators
H = rand_hermitian(dim);
D = rand_hermitian(dim);
[dH, A] = markov.ops(H, D);

X = A{1}; % dH(1) == 0
for k=2:length(dH)
  X = X +A{k} +A{k}'; % A(-omega) == A'(omega)
end
assert_o(dH(1), 0, tol);
assert_o(norm(X-D), 0, tol); % Lindblad ops should sum to D
