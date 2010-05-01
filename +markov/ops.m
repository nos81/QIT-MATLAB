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

% Ville Bergholm 2009-2010


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
p = length(deltaE);


if (~iscell(D))
  D = {D}; % D needs to be a cell array, even if it has just one element
end
n_D = length(D); % number of bath coupling ops

% combine degenerate deltaE, build jump ops
% k -> ind(k) -> i,j
s = 1;
[r,c] = ind2sub([m m], ind(1));
dH(s) = deltaE(1);
  
for op=1:n_D
  A{op,s} = P{c} * D{op} * P{r};
end
  
for k = 2:p
  [r,c] = ind2sub([m m], ind(k));

  if (abs(deltaE(k) - deltaE(k-1)) > qit.tol)
    % new omega value, new jump op
    s = s+1;
    dH(s) = deltaE(k);
    for op=1:n_D
      A{op,s} = P{c} * D{op} * P{r};
    end
  else
    % extend current op
    for op=1:n_D
      A{op,s} = A{op,s} +P{c} * D{op} * P{r};
    end
  end
end
