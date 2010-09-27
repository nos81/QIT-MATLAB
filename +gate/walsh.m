function U = walsh(n)
% WALSH  Walsh-Hadamard gate.
%  U = walsh(n)
%
%  Returns the Walsh-Hadamard gate for n qubits.

% Ville Bergholm 2009-2010


H = [1 1; 1 -1]/sqrt(2);
U = 1;
for k = 1:n
  U = kron(U, H);
end

dim = qubits(n);
U = lmap(U, {dim, dim});
