function [f] = fermion_ladder(n)
% FERMION_LADDER  Fermionic ladder operators.
%  f = fermion_ladder(n)
%
%  Returns a cell vector of fermionic annihilation operators for a
%  system of n fermionic modes in the second quantization.
%
%  The annihilation operators are built using the Jordan-Wigner
%  transformation for a chain of n qubits, where the state of each
%  qubit denotes the occupation number of the corresponding mode.
%
%  First define annihilation and number operators for a lone fermion mode:
%   s := (sx + i*sy)/2   = [0 1; 0 0],
%   n := s'*s = (I-sz)/2 = [0 0; 0 1].
%
%   s|0> = 0, s|1> = |0>, n|k> = k|k>
%
%  Then define a phase operator to keep track of sign changes when
%  permuting the order of the operators:
%   \phi_k := \sum_{j=1}^{k-1} n_j.
%
%  Now, the fermionic annihilation operators for the n-mode system are given by
%   f_k := (-1)^{\phi_k} s_k.
%
%  These operators fulfill the required anticommutation relations:
%   {f_k, f_j}  = 0,
%   {f_k, f_j'} = I \delta_{kj},
%   f_k' * f_k  = n_k.

% Ville Bergholm 2009-2010


d = 2^n;

% number and phase operators (diagonal, so we store them as such)
temp = zeros(d, 1);
phi{1} = temp;
for k=1:n-1
  num = mkron(ones(2^(k-1), 1), [0; 1], ones(2^(n-k), 1)); % number operator n_k as a diagonal
  temp = temp + num; % sum of number ops up to n_k, diagonal
  phi{k+1} = temp;
end

s = sparse([0 1; 0 0]);

% fermionic annihilation operators
for k=1:n
  f{k} = spdiags((-1).^phi{k}, 0, d, d) * mkron(speye(2^(k-1)), s, speye(2^(n-k)));
end
