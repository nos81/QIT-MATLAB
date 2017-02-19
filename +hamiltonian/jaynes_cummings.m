function [H, dim, H_int] = jaynes_cummings(om_atom, om_cavity, J, m, use_RWA)
% JAYNES_CUMMINGS  Jaynes-Cummings model, two-level atom coupled to a cavity.
%  [H, dim] = jaynes_cummings(om_atom, om_cavity, J, m, use_RWA)
%
%  Returns the Hamiltonian H and the dimension vector dim for an
%  implementation of the Jaynes-Cummings model, describing n two-level atoms coupled
%  to c harmonic oscillators (e.g. individual EM field modes in an optical cavity),
%  where n == length(om_atom), c == length(om_cavity), and size(J) == [c, n].
%
%  H/\hbar = -\sum_k \frac{\omega_a_k}{2} \sigma_z_k +\sum_j \omega_c_j a_j^\dagger a_j +\sum_{jk} \frac{J_{jk}}{2} \sigma_x_k (a_j+a_j^\dagger)
%
%  The Hamiltonian H can also be normalized with e.g. \omega_c_1, and thus be made dimensionless.
%
%  The order of the subsystems is [cavity_1, ..., cavity_c, atom_1, ..., atom_n].
%  The dimension of the Hilbert space of each bosonic cavity mode (infinite in principle) is truncated to m.
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
% number of cavity modes
c = length(om_cavity);

if size(J) ~= [c, n]
    error('The coupling matrix J must be of size [length(om_cavity), length(om_atom)].');
end

% dimension vector
dim = [m*ones(1,c), 2*ones(1,n)];

% operators
a = boson_ladder(m);
x = a+a';
sp = 0.5*(qit.sx -1i*qit.sy); % qubit raising operator

atom = cell(1, n);
cavity = cell(1, c)
coupling = cell(c, n);
% loop over cavity modes
for j=1:c
    cavity{j} = {om_cavity(j) * a' * a, j}
end
% loop over atoms
for k=1:n
    atom{k} = {-0.5*om_atom(k) * qit.sz, k+c}; % atomic Hamiltonian

    % loop over cavity modes
    for j=1:c
        % atom-cavity coupling
        if use_RWA
            % rotating wave approximation, discard counter-rotating terms
            coupling{j,k}   = {a,  j; 0.5*J(j,k) * sp,  k+c};
            coupling{j,n+k} = {a', j; 0.5*J(j,k) * sp', k+c};
        else
            coupling{j,k} = {x, j; 0.5*J(j,k) * qit.sx, k+c};
        end
    end
end
Ha = op_list(atom, dim);
Hc = op_list(cavity, dim);
H_int = op_list(coupling, dim);

if nargout == 3
    % return interaction terms separately
    H = Hc +Ha;
else
    H = Hc +Ha +H_int;
end
