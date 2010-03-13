function [H, dim] = hubbard(n, U_t, mu_t)
% HUBBARD  Hubbard model, fermions in a 1D lattice.
%  [H, dim] = hubbard(n, U_t, mu_t)
%
%  Returns the Hamiltonian H and the dimension vector dim for an
%  implementation of the Hubbard model with n lattice sites.
%
%  The model consists of spin-1/2 fermions confined in a
%  one-dimensional lattice. The fermions interact with other fermions
%  in the same site, as well as with an external chemical potential.
%
%  H = \sum_k (-\sum_s f_{k,s}' f_{k+1,s} +h.c.) +U/t n_{k,1} n_{k,2} -\mu/t (n_{k,1}+n_{k,2})
%
%  The Hamiltonian has been normalized by the fermion hopping
%  constant t. The other parameters are U_t == U/t and mu_t = \mu/t. 

% Ville Bergholm 2010


dim = 2*ones(1, 2*n); % n sites, two fermionic modes per site

% fermion annihilation ops (same-site states next to each other, f{site,spin})
f = fermion_ladder([2, n]).'; 

if (nargin < 3)
  mu_t = 0;
  if (nargin < 2)
    U_t = 1;
  end
end

H = sparse(0);

for k=1:n-1
  % fermions hopping
  for s=1:2
    H = H -(f{k,s}'*f{k+1,s} +f{k+1,s}'*f{k,s});
  end
end

for k=1:n
  % number operators for this site
  n1 = f{k,1}'*f{k,1};
  n2 = f{k,2}'*f{k,2};
  
  % on-site interaction
  H = H +U_t * n1 * n2;
  
  % chemical potential
  H = H -mu_t * (n1 + n2);
end
