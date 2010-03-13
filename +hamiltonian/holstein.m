function [H, dim] = holstein(n, m, omega_t, g)
% HOLSTEIN  Holstein model, electrons in a 1D lattice coupled to phonons.
%  [H, dim] = holstein(n, m, omega_t, g)
%
%  Returns the Hamiltonian H and the dimension vector dim for an
%  implementation of the Holstein model with n lattice sites.
%
%  The model consists of spinless electrons confined in a
%  one-dimensional lattice, coupled to phonon modes represented
%  by a harmonic oscillator at each site. The dimensions of phonon
%  Hilbert spaces (infinite in principle) are truncated to m.
%
%  The order of the subsystems is [e1, ..., en, p1, ..., pn].
%
%  H = \sum_k (-c_k' c_{k+1} +h.c.) +\omega/t b_k' b_k -g \omega/t (b_k + b_k') c_k' c_k
%
%  The Hamiltonian has been normalized by the electron hopping
%  constant t. The other parameters are omega_t == \omega/t and g. 

% Ville Bergholm 2010


% Hilbert space: electrons first, then phonons
dim = [2^n, m*ones(1, n)] % j-w clumps fermion dims together

b = ho.ladder(m); % phonon annihilation
q = b+b';  % phonon position
nb = b'*b; % phonon number operator

c = fermion_ladder(n); % electron annihilation ops

if (nargin < 4)
  omega_t = 1;
  g = 1;
end

H = sparse(0);

for k=1:n
  % phonon harmonic oscillators
  H = H +omega_t * op_list({{nb, 1+k}}, dim);

  % electron-phonon interaction
  H = H -g * omega_t * op_list({{c{k}'*c{k},1; q,1+k}}, dim);
end

T = sparse(0);
for k=1:n-1
  % electrons hopping
  T = T -(c{k}'*c{k+1} +c{k+1}'*c{k});
end
H = H+kron(T, speye(prod(dim(2:end))));

% actual dimensions
dim = [2*ones(1, n), m*ones(1, n)];
