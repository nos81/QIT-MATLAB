function [H, dim, H_int] = jaynes_cummings(om_atom, Om, m, use_RWA)
% JAYNES_CUMMINGS  Jaynes-Cummings model, one or more two-level atoms coupled to a single-mode cavity.
%  [H, dim] = jaynes_cummings(om_atom, Om, m, use_RWA)
%
%  Returns the Hamiltonian H and the dimension vector dim for an
%  implementation of the Jaynes-Cummings model, describing n two-level atoms coupled
%  to a harmonic oscillator (e.g. a single EM field mode in an optical cavity),
%  where n = length(om_atom) = length(Om).
%
%  H/\hbar = -\sum_k \frac{\omega_a_k}{2} \sigma_z_k +\omega_c a^\dagger a +\sum_k \frac{\Omega_k}{2} \sigma_x_k (a+a^\dagger)
%
%  The returned Hamiltonian H has been additionally normalized with \omega_c,
%  and is thus dimensionless. om_atom(k) = \omega_a_k / \omega_c,  Om(k) = \Omega_k / \omega_c.
%
%  The order of the subsystems is [cavity, atom_1, ..., atom_n].
%  The dimension of the Hilbert space of the bosonic cavity mode (infinite in principle) is truncated to m.
%  If use_RWA is true, the Rotating Wave Approximation is applied to the Hamiltonian,
%  and the counter-rotating interaction terms are discarded.

% Ville Bergholm 2014-2016


if nargin < 3
    m = 10;
end
if nargin < 4
    use_RWA = 0;
end

global qit

% number of atoms
n = length(om_atom);
if length(Om) ~= n
    error('The coupling vector Om must be of the same length as the atom splitting vector om_atom.');
end

% dimension vector
dim = [m, 2*ones(1,n)];

% operators
a = boson_ladder(m);
x = a+a';
sp = 0.5*(qit.sx -1i*qit.sy); % qubit raising operator

% cavity
Hc = op_list({{a' * a, 1}}, dim);

atom = cell(1, n);
coupling = cell(1, n);
% loop over atoms
for k=1:n
    atom{k} = {-0.5*om_atom(k) * qit.sz, k+1}; % atomic Hamiltonian
    % atom-cavity coupling
    if use_RWA
        % rotating wave approximation, discard counter-rotating terms
        coupling{k}   = {a,  1; 0.5*Om(k) * sp,  k+1};
        coupling{n+k} = {a', 1; 0.5*Om(k) * sp', k+1};
    else
        coupling{k} = {x, 1; 0.5*Om(k) * qit.sx, k+1};
    end
end
Ha = op_list(atom, dim);
H_int = op_list(coupling, dim);

if nargout == 3
    % return interaction terms separately
    H = Hc +Ha;
else
    H = Hc +Ha +H_int;
end
