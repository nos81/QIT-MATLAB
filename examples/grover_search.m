function [p] = grover_search(n)
% GROVER_SEARCH  Grover search algorithm demo.
%
%  p = grover_search(n)
%
%  Simulate the Grover search algorithm formulated using amplitude amplification
%  in a system of n qubits.

%! L.K. Grover, "Quantum Computers Can Search Rapidly by Using Almost Any Transformation", PRL 80, 4329 (1998). doi:10.1103/PhysRevLett.80.4329.   
% Ville Bergholm 2009


fprintf('\n\n=== Grover search algorithm ===\n')

global qit;

A = gate.walsh(n); % Walsh-Hadamard gate for generating uniform superpositions
N = 2^n; % number of states

sol = floor(rand() * N) + 1; % 1..N
reps = floor(pi/(4*asin(sqrt(1/N))));

fprintf('Using %d qubits\n', n)
fprintf('Solution: %d\n', sol)
fprintf('Probability maximized by %d iterations\n', reps)


% black box oracle capable of recognizing the correct answer (given as the diagonal)
% TODO an oracle that actually checks the solutions by computing (using ancillas?)
U_oracle = ones(N, 1);
U_oracle(sol) = -1;



s = state(0, 2*ones(n, 1));

% initial superposition
s = u_propagate(s, A);

% Grover iteration
for k=1:reps
  % oracle phase flip
  s.data = -U_oracle .* s.data;

  % inversion about the mean
  s = u_propagate(s, A');
  temp = s.data; % FIXME annoying
  temp(1) = -temp(1); % phase flip the zero state
  s.data = temp;
  s = u_propagate(s, A);
end

p = prob(s);
