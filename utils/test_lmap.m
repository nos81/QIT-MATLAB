% Test script for the lmap class.
% Ville Bergholm 2009-2010

disp('Testing the lmap class.')

tol = qit.tol;

% reordering subsystems
dim = [2 5 3];
A = lmap(rand(dim(1)));
B = lmap(rand(dim(2)));
C = lmap(rand(dim(3)));
T1 = tensor(A, B, C);

p = [3 1 2];
T2 = reorder(T1, {p, p});
assert_o(norm(tensor(C, A, B) - T2), 0, tol);

p = [2 1 3];
T2 = reorder(T1, {p, p});
assert_o(norm(tensor(B, A, C) - T2), 0, tol);

p = [3 2 1];
T2 = reorder(T1, {p, p});
assert_o(norm(tensor(C, B, A) - T2), 0, tol);





% Test script for gates.
% Ville Bergholm 2010

%U = gate.two(B, t, dim);
