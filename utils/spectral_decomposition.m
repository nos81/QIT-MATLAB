function [a, P] = spectral_decomposition(A, do_sort, tol)
% SPECTRAL_DECOMPOSITION  Spectral decomposition of a Hermitian matrix.
%  [a, P] = spectral_decomposition(A[, sort, tol])
%
%  Returns the unique eigenvalues a and the corresponding projectors P
%  for the Hermitian matrix A, such that  A = \sum_k a_k P_k.
%  If sort is true, the eigenvalues are sorted in ascending order.
%  If tol is given, it determines the tolerance for grouping
%  eigenvalues close together.

% Ville Bergholm 2010-2017


global qit
if nargin < 3
    % tolerance for grouping eigenvalues together
    tol = qit.tol;
end

[v,d] = eig(A);
d = real(diag(d)); % A is assumed Hermitian

if nargin >= 2 && do_sort
    [d, ind] = sort(d); % ascending order
    v = v(:, ind);
end

% combine projectors for degenerate eigenvalues
s = 1;
a(1,1) = d(1);
P{1} = v(:,1)*v(:,1)';

for k = 2:length(d)
  if abs(d(k) - d(k-1)) > tol
    % new eigenvalue, new projector
    s = s+1;
    a(s,1) = d(k);
    P{s} = v(:,k)*v(:,k)';
  else
    % extend current P
    P{s} = P{s} + v(:,k)*v(:,k)';
  end
end
