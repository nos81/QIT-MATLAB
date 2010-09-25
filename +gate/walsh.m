function U = walsh(n)
% WALSH  Walsh-Hadamard gate.
%  U = walsh(n)
%
%  Returns the Walsh-Hadamard gate for n qubits.

% Ville Bergholm 2009

H = [1 1; 1 -1]/sqrt(2);
U = 1;
for k = 1:n
  U = kron(U, H);
end

dim = 2*ones(1, n);
U = lmap(U, {dim, dim});
