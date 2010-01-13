function [dH, A] = ops(H, D)
% LINDBLAD/OPS  Lindblad operators for a Born-Markov master equation.
%  [dH, A] = ops(H, D)
%
%  Builds the Lindblad operators for a Hamiltonian operator H and
%  a (hermitian) interaction operator D.
%
%  Returns dH, a vector of the unique nonnegative differences between
%  eigenvalues of H, and A, an array of the corresponding Lindblad operators.
%  length(dH) == length(A)

% Ville Bergholm 2009


tol = 1e-12;

[E, P] = spectral_decomposition(H);
m = length(E); % unique eigenvalues

%X = 0;
%for k=1:m
%  X = X + E(k)*P{k};
%end
%assert(norm(X-H), 0, tol);


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

% combine degenerate deltaE, build Lindblad ops
% k -> ind(k) -> i,j
s = 1;
[r,c] = ind2sub([m m], ind(1));
dH(s) = deltaE(1);
  
for op=1:n_D
  A{op,s} = P{c} * D{op} * P{r};
end
  
for k = 2:p
  [r,c] = ind2sub([m m], ind(k));

  if (abs(deltaE(k) - deltaE(k-1)) > tol)
    % new omega value, new Lindblad op
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
