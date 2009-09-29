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

n = length(H); % total eigenvalues

[v,d] = eig(H);
d = diag(d);

% projectors to energy subspaces
s = 1;
E(1) = d(1);
P{1} = v(:,1)*v(:,1)';

% get rid of degenerate energies
for k = 2:n
  if (abs(d(k) - d(k-1)) > tol)
    % new energy value, new projector
    s = s+1;
    E(s) = d(k);
    P{s} = v(:,k)*v(:,k)';
  else
    % extend current projector
    P{s} = P{s} + v(:,k)*v(:,k)';
  end
end

m = length(E); % unique eigenvalues
% now H = \sum_k E_k P_k

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


n_D = length(D); % number of bath coupling ops
for op=1:n_D

  % combine degenerate deltaE, build Lindblad ops
  % k -> ind(k) -> i,j
  s = 1;
  [r,c] = ind2sub([m m], ind(1));
  dH(s) = deltaE(1);
  A{op,s} = P{c} * D{op} * P{r};

  for k = 2:p
    [r,c] = ind2sub([m m], ind(k));

    if (abs(deltaE(k) - deltaE(k-1)) > tol)
      % new omega value, new Lindblad op
      s = s+1;
      dH(s) = deltaE(k);
      A{op,s} = P{c} * D{op} * P{r};
    else
      % extend current op
      A{op,s} = A{op,s} +P{c} * D{op} * P{r};
    end
  end
end
