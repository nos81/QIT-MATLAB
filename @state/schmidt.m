function [lambda, v, u] = schmidt(s, sys)
% SCHMIDT  Schmidt decomposition.
%  lambda = schmidt(s, sys)
%  [lambda, u, v] = schmidt(s, sys)
%
%  Calculates the Schmidt decomposition of the (pure) state s.
%  Subsystems listed in vector sys constitute part A, the rest forming part B.
%  Vector lambda will contain the Schmidt coefficients.
%
%  If required, matrices u and v will contain the corresponding orthonormal
%  Schmidt bases for A and B, respectively, as column vectors, i.e.
%  \ket{k}_A = u(:,k), \ket{k}_B = v(:,k).
%  The state s is then given by \sum_k \lambda_k \ket{k}_A \otimes \ket{k}_B

%! M.A. Nielsen, I.L. Chuang, "Quantum Computation and Quantum Information" (2000), chapter 2.5.
% Ville Bergholm 2009-2010


% number of systems
dim = dims(s);
n = length(dim);

if (nargin < 2)
  if (n == 2)
    % reasonable choice
    sys = 1;
  else
    error('Requires a state and a vector of subsystems.')
  end
end

if ~is_ket(s)
  % state operator
  s = to_ket(s);
  %error('Schmidt decomposition is only defined for pure states.')
end

% complement of sys, dimensions of the partitions
sys = clean_selection(s, sys);
compl = setdiff(1:n, sys);
d1 = prod(dim(sys));
d2 = prod(dim(compl));
perm = [sys, compl];

if all(perm == 1:n)
  % do nothing
else
  % reorder the system according to the partitioning
  s = reorder(s, perm);
end

% order the coefficients into a matrix, take an svd (reshape wants little-endian dims)
if (nargout == 1)
  lambda = svd(reshape(s.data, [d2 d1]), 'econ');
else
  [u, lambda, v] = svd(reshape(s.data, [d2 d1]), 'econ');
  v = conj(v); % note the definition of v in svd
  lambda = diag(lambda);
end
