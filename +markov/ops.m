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

% Ville Bergholm 2009-2015


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

s = 0; % index of the currently built jump op
current_dH = Inf;

% combine degenerate deltaE, build jump ops
for k = 1:length(deltaE)
  diff = abs(deltaE(k) -current_dH);
  if diff > qit.tol
    % new omega value, new jump op
    % otherwise just extend current op
    s = s+1;
    dH(s) = deltaE(k);
    current_dH = dH(s);
    A(:,s) = num2cell(zeros(n_D, 1));
  end
  % given k, find corresponding row and column in the deltaE matrix
  [r, c] = ind2sub([m, m], ind(k));
  % add projection of D corresponding to the index k to the jump op number s
  for op=1:n_D
      A{op, s} = A{op, s} +P{c} * D{op} * P{r};
  end
end

% eliminate zero As and corresponding dHs
temp = zeros(size(A));
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
% TODO justify the numerical tolerance used
for k=2:length(dH)
    if abs(dH(k)-dH(k-1)) < 1e-3
        fprintf('Warning: Small difference between dH(%d) and dH(%d) may break the RWA.\n', k-1, k);
    end
end
end
