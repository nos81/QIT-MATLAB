function adiabatic_qc(n, n_clauses)
% ADIABATIC_QC  Adiabatic quantum computing example.
%  adiabatic_qc(n_bits, n_clauses)
%
%  This example solves random 3-SAT problems using the
%  adiabatic quantum algorithm due to Farhi et al.
%
%  Note that this is incredibly inefficient because we first essentially
%  solve the NP-complete problem using an exhaustive search when computing
%  the problem Hamiltonian H1, and then simulate an adiabatic quantum computer
%  solving the same problem using the quantum algorithm.

%! E. Farhi et al., "Quantum Computation by Adiabatic Evolution", arXiv.org:quant-ph/0001106.
% Ville Bergholm 2009


fprintf('\n\n=== Solving 3-SAT using adiabatic qc ===\n')

if (nargin < 2)
  if (nargin < 1)
    n = 6;
  end
  n_clauses = n/2;
end

if (n < 3)
  n = 3;
end

% generate clauses, count how many times each bit is referenced in the clauses
clauses = zeros(n_clauses, 3);
refs = zeros(n, 1);
for k=1:n_clauses
  bits = 1:n;
  for j=1:3
    choice = randi(length(bits));
    temp = bits(choice);
    bits(choice) = [];
    clauses(k, j) = temp;
    refs(temp) = refs(temp) + 1;
  end
  clauses(k, :) = sort(clauses(k, :));
end

disp('Find a bit string such that all the following "exactly one of bits {a, b, c} is 1" clauses are satisfied:')
clauses

% build initial Hamiltonian
global qit;
xb = 0.5*(qit.I - qit.sx);
H0 = 0;
for k=1:n
  H0 = H0 + refs(k)*mkron(eye(2^(k-1)), xb, eye(2^(n-k)));
end

% initial state (ground state of H0)
s = state(ones(2^n, 1)/sqrt(2^n), 2*ones(1, n)); % n qubits, uniform superposition

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
  % h_c = (b1 ^ b2* ^ b3*) v (b1* ^ b2 ^ b3*) v (b1* ^ b2* ^ b3)
  b1 = z_op{clauses(k, 1)};
  b2 = z_op{clauses(k, 2)};
  b3 = z_op{clauses(k, 3)};
  s1 = b1.*not(b2).*not(b3);
  s2 = b2.*not(b3).*not(b1);
  s3 = b3.*not(b1).*not(b2);
  H1 = H1 +not(s1+s2+s3 -s1.*s2 -s2.*s3 -s3.*s1 +s1.*s2.*s3);
end
H1_full = diag(H1); % into a full matrix

% adiabatic simulation
t = 50;
steps = t*10;
res = adiabatic_propagate(s, H0, H1_full, t, steps);

% plots
% final state probabilities
figure;
plots.tomography(res{end});
title('Final state');

H1

disp('Measured result:')
[p, dummy, res] = measure(res{end});
display(res)
if (H1(find(res.data)) == 0)
  disp('Which is a valid solution!')
else
  disp('Which is not a solution!')
end


end

function a = randi(n)
  a = ceil(rand*n);
end
