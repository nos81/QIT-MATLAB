function [H1, clauses] = adiabatic_qc_3sat(n, n_clauses, clauses)
% ADIABATIC_QC_3SAT  Solving 3-SAT using adiabatic quantum computing.
%  [H1, clauses] = adiabatic_qc_3sat(n_bits, n_clauses)
%
%  This example solves random 3-SAT problems by simulating the
%  adiabatic quantum algorithm of Farhi et al.
%
%  Note that this is incredibly inefficient because we first essentially
%  solve the NP-complete problem classically using an exhaustive
%  search when computing the problem Hamiltonian H1, and then
%  simulate an adiabatic quantum computer solving the same problem
%  using the quantum algorithm.

%! E. Farhi et al., "Quantum Computation by Adiabatic Evolution", arXiv.org:quant-ph/0001106.
% Ville Bergholm 2009-2010


fprintf('\n\n=== Solving 3-SAT using adiabatic qc ===\n\n')

if (nargin < 2)
  if (nargin < 1)
    n = 6;
  end
  n_clauses = 5*n;
end

if (n < 3)
  n = 3;
end

fprintf('%d bits, %d clauses\n', n, n_clauses);

if (nargin < 3)
% generate clauses
clauses = zeros(n_clauses, 3);
for k=1:n_clauses
  bits = 1:n;
  for j=1:3
    choice = randi(length(bits));
    temp = bits(choice);
    bits(choice) = [];
    clauses(k, j) = ((-1)^(rand < 0.5))*temp; % negate if bit should be inverted
  end
  [~, I] = sort(abs(clauses(k, :)));
  clauses(k, :) = clauses(k, I);
end
end

disp('Find a bit string such that all the following "(b_i or b_j or b_k)" clauses are satisfied.')
disp('Minus sign means the bit is negated.')
clauses

% cache some stuff (all the matrices in this example are diagonal, so)
zb  = [0 1]; % 0.5*(I - sz);
z_op = [];
for k=1:n
  z_op{k} = mkron(ones(1, 2^(k-1)), zb, ones(1, 2^(n-k)));
end

% Encode a 3-SAT problem into the final Hamiltonian.
% The clauses we use here correspond to the exact cover problem.
H1 = 0;
for k=1:n_clauses
  % h_c = (b1(*) v b2(*) v b3(*))
  for j=1:3
    b{j} = z_op{abs(clauses(k, j))};
    if (clauses(k, j) < 0)
      b{j} = not(b{j});
    end
  end
  H1 = H1 +not(b{1}+b{2}+b{3}); % the not makes this a proper OR
end

% build initial Hamiltonian
global qit;
xb = 0.5*(qit.I - qit.sx);
H0 = 0;
for k=1:n
    H0 = H0 +mkron(eye(2^(k-1)), xb, eye(2^(n-k)));
end

% initial state (ground state of H0)
s0 = state(ones(2^n, 1)/sqrt(2^n), 2*ones(1, n)); % n qubits, uniform superposition


% adiabatic simulation
adiabatic_qc(H0, H1, s0);
end

function a = randi(n)
  a = ceil(rand*n);
end
