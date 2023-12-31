function [p] = grover_search(n)
% GROVER_SEARCH  Grover search algorithm demo.
%
%  p = grover_search(n)
%
%  Simulate the Grover search algorithm formulated using amplitude amplification
%  in a system of n qubits.

%! L.K. Grover, "Quantum Computers Can Search Rapidly by Using Almost Any Transformation", PRL 80, 4329 (1998). doi:10.1103/PhysRevLett.80.4329.   
% Ville Bergholm 2009-2010


fprintf('\n\n=== Grover search algorithm ===\n')

global qit;

A = gate.walsh(n); % Walsh-Hadamard gate for generating uniform superpositions
N = 2^n; % number of states

sol = randi(N);
reps = floor(pi/(4*asin(sqrt(1/N))));

fprintf('Using %d qubits.\n', n)
fprintf('Probability maximized by %d iterations.\n', reps)
fprintf('Correct solution: %d\n', sol)


% black box oracle capable of recognizing the correct answer (given as the diagonal)
% TODO an oracle that actually checks the solutions by computing (using ancillas?)
U_oracle = ones(N, 1);
U_oracle(sol) = -1;

U_zeroflip = ones(N, 1);
U_zeroflip(1) = -1;

s = state(0, qubits(n));

% initial superposition
s = u_propagate(s, A);

% Grover iteration
for k=1:reps
  % oracle phase flip
  s.data = -U_oracle .* s.data;

  % inversion about the mean
  s = u_propagate(s, A');
  s.data = U_zeroflip .* s.data;
  s = u_propagate(s, A);
end

[dummy, res] = measure(s);
fprintf('\nMeasured %d.\n', res);
