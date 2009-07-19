function U = walsh(n)
% walsh  Walsh-Hadamard gate.
%  U = walsh(n)
%
%  Returns the Walsh-Hadamard matrix for n qubits.

% Ville Bergholm 2009

H = [1 1; 1 -1]/sqrt(2);
U = 1;
for k = 1:n
  U = kron(U, H);
end
