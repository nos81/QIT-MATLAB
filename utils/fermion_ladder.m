function [f] = fermion_ladder(grouping)
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
%   f_k := (-1)^{\phi_k} s_k = (\bigotimes_{j=1}^{k-1} sz_j) s_k.
%
%  These operators fulfill the required anticommutation relations:
%   {f_k, f_j}  = 0,
%   {f_k, f_j'} = I \delta_{kj},
%   f_k' * f_k  = n_k.

% Ville Bergholm 2009-2016


if isscalar(grouping)
    grouping = [1, grouping];
end

n = prod(grouping);

s  = sparse([0, 1; 0, 0]);  % single annihilation op
sz = sparse([1, 0; 0,-1]);

% empty cell array for the annihilation operators
f = cell(n,1);

% Jordan-Wigner transform
temp = 1;
for k=1:n
    f{k} = mkron(temp, s, speye(2^(n-k)));
    temp = kron(temp, sz);
end

f = reshape(f, grouping);
