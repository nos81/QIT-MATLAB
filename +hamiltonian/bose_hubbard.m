function [H, dim] = bose_hubbard(n, m, U_t, mu_t)
% BOSE_HUBBARD  Bose-Hubbard model, bosons in a 1D lattice.
%  [H, dim] = bose_hubbard(n, m, U_t, mu_t)
%
%  Returns the Hamiltonian H and the dimension vector dim for an
%  implementation of the Bose-Hubbard model with n lattice sites.
%
%  The model consists of spinless bosons confined in a
%  one-dimensional lattice. The bosons interact with other bosons
%  in the same site, as well as with an external chemical potential.
%  The dimensions of boson Hilbert spaces (infinite in principle)
%  are truncated to m.
%
%  H = \sum_k -(b_k' b_{k+1} +h.c.) +U/(2t) n_k (n_k-1) -\mu/t n_k
%
%  The Hamiltonian has been normalized by the boson hopping
%  constant t. The other parameters are U_t == U/t and mu_t = \mu/t. 

% Ville Bergholm 2010


dim = m*ones(1, n);

b = ho.ladder(m); % boson annihilation
nb = b'*b; % boson number operator

if (nargin < 4)
  mu_t = 0;
  if (nargin < 3)
    U_t = 1;
  end
end

I = speye(m);
A = U_t/2 * nb * (nb-I); % on-site interaction
B = -mu_t * nb; % chemical potential

H = sparse(0);

for k=1:n
  H = H +op_list({{A+B, k}}, dim);
end

for k=1:n-1
  % bosons hopping
  H = H -op_list({{b', k; b, k+1}, {b, k; b', k+1}}, dim);
end
