function [dH, A] = ops(H, D)
% OPS  Jump operators for a Born-Markov master equation.
%  [dH, A] = ops(H, D)
%
%  Builds the jump operators for a Hamiltonian operator H and
%  a (hermitian) interaction operator D.
%
%  Returns dH, a vector of the sorted unique nonnegative differences between
%  eigenvalues of H, and A, a cell array of the corresponding jump operators.
%  size(A) == [length(D), length(dH)]

% Ville Bergholm 2009-2010, 2015


global qit;

[E, P] = spectral_decomposition(full(H));
m = length(E); % unique eigenvalues

% energy difference matrix is antisymmetric, so we really only need the lower triangle
deltaE = kron(E.', ones(1,m)) -kron(E, ones(m,1)); % deltaE(i,j) = E(i)-E(j)

% sort is a stable sorting algorithm
[deltaE, ind] = sort(deltaE(:));

% index of first lower triangle element
s = m*(m-1)/2 + 1;
%assert(ind(s), 1)
deltaE = deltaE(s:end).';
ind = ind(s:end).';

if ~iscell(D)
  D = {D}; % D needs to be a cell array, even if it has just one element
end
n_D = length(D); % number of bath coupling ops

% jump ops
A = cell(n_D, 0);

s = 1; % index of the currently built jump op
% start the first jump op
add_projection(1, s);
dH(s) = deltaE(1);

% combine degenerate deltaE, build jump ops
for k = 2:length(deltaE)
  diff = abs(deltaE(k) - dH(s));
  if diff > qit.tol
    % new omega value, new jump op
    s = s+1;
    dH(s) = deltaE(k);
    add_projection(k, s);
  else
    % extend current op
    add_projection(k, s);
  end
end

    function add_projection(X_k, X_s)
    % Adds a projection of D corresponding to the index X_k to the
    % jump op number X_s.
    % NOTE shares scope with parent function! We use the X_ prefix
    % for local variables to avoid accidentally overwriting parent
    % function stuff.
    [X_r, X_c] = ind2sub([m, m], ind(X_k));

    for X_op=1:n_D
        A{X_op, X_s} = P{X_c} * D{X_op} * P{X_r};
    end
    end

% eliminate zero As and corresponding dHs
temp = [];
for k=1:length(dH)
    for op=1:n_D
        temp(op, k) = norm(A{op, k}) < qit.tol;
    end
end
% columns in which every A vanishes
temp = all(temp, 1);
% remove them
A(:,temp) = [];
dH(temp) = [];

% Are some of the remaining dH differences too low for RWA to hold properly?
for k=2:length(dH)
    if abs(dH(k)-dH(k-1)) < 1e-3
        fprintf('Warning: Small difference between dH(%d) and dH(%d) may break the RWA.\n', k-1, k);
    end
end
end
