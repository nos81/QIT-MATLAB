function [H, dim] = heisenberg(dim, C, J, B)
% HEISENBERG  Heisenberg spin network model
%  [H, dim] = heisenberg(dim, C, J, B)
%
%  Returns the Hamiltonian H for the Heisenberg model, describing a
%  network of n interacting spins in an external magnetic field.
%
%  dim is a dimension vector of the spins, i.e. dim == [2 2 2]
%  would be a system of three spin-1/2's.
%
%  C is the n \times n connection matrix of the spin network, where C(i, j)
%  is the coupling strength between spins i and j. Only the upper triangle is used.
%
%  J defines the form of the spin-spin interaction. It is either a
%  3-vector (homogeneous couplings) or a function J(i, j) returning
%  a 3-vector for site-dependent interactions.
%  Element k of the vector is the coefficient of the Hamiltonian term S_k^{(i)} S_k^{(j)},
%  where S_k^{(i)} is the k-component of the angular momentum of spin i.
%
%  B defines the effective magnetic field the spins locally couple to. It's either
%  a 3-vector (homogeneous field) or a function B(i) that returns a 3-vector for
%  site-dependent field.
%
%  Examples:
%    C = diag(ones(1, n-1), 1)  linear n-spin chain
%    J = [2 2 2]                isotropic Heisenberg coupling
%    J = [2 2 0]                XX+YY coupling
%    J = [0 0 2]                Ising ZZ coupling
%    B = [0 0 1]                homogeneous Z-aligned field
%
%  H =  \sum_{\langle i,j \rangle} \sum_{k = x,y,z} J(i,j)[k] S_k^{(i)} S_k^{(j)}
%      +\sum_i \vec{B}(i) \cdot \vec{S}^{(i)})

% Ville Bergholm 2009-2014


if nargin < 1
    error('Usage: heisenberg(dim, C, J, B)')
end
n = length(dim);

if nargin < 2
    % linear open chain
    C = 'open';
end
if nargin < 3
    % Ising ZZ
    J = [0, 0, 2];
end
if nargin < 4
    % Z field
    B = [0, 0, 1];
end

% simple way of defining some common cases
if ischar(C)
    switch C
      case 'open'
        C = diag(ones(1, n-1), 1);
      case 'closed'
        C = diag(ones(1, n-1), 1);
        C(1, end) = 1;
      otherwise
        error('Unknown topology.')
    end
end

if isa(J, 'function_handle')
    Jfunc = J;
else
    if isvector(J) && length(J) == 3
        Jfunc = @(i, j) C(i, j) * J;
    else
        error('J must be either a 3x1 vector or a function handle.')
    end
end

if isa(B, 'function_handle')
    Bfunc = B;
else
    if isvector(B) && length(B) == 3
        Bfunc = @(i) B;
    else
        error('B must be either a 3x1 vector or a function handle.')
    end
end

% local magnetic field terms
temp = {};
for i = 1:n
    A = angular_momentum(dim(i));  % spin ops
    temp = cat(2, temp, {{cdot(Bfunc(i), A), i}});
end
H = op_list(temp, dim);

% spin-spin couplings: loop over nonzero entries of C
% only use the upper triangle
[a, b] = find(triu(C));
for m = 1:length(a)
    i = a(m);
    j = b(m);
    % spin ops for sites i and j
    Si = angular_momentum(dim(i));
    Sj = angular_momentum(dim(j));
    temp = {};
    % coupling between sites a and b
    c = Jfunc(i, j);
    for k = 1:3
        temp = cat(2, temp, {{c(k) * Si{k}, i; Sj{k}, j}});
        H = H +op_list(temp, dim);
    end
end
end

function res = cdot(v, A)
% Real dot product of a vector with a cell vector of operators.

res = 0;
for k = 1:length(v)
    res = res +v(k) * A{k};
end
end
