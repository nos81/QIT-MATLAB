function [a, P] = spectral_decomposition(A, do_sort)
% SPECTRAL_DECOMPOSITION  Spectral decomposition of a Hermitian matrix.
%  [a, P] = spectral_decomposition(A[, sort])
%
%  Returns the unique eigenvalues a and the corresponding projectors P
%  for the Hermitian matrix A, such that  A = \sum_k a_k P_k.
%  If sort is true, the eigenvalues are sorted in ascending order.

% Ville Bergholm 2010-2013


global qit

[v,d] = eig(A);
d = real(diag(d)); % A is assumed Hermitian

if nargin >= 2 && do_sort
    [d, ind] = sort(d); % ascending order
    v = v(:, ind);
end

% combine projectors for degenerate eigenvalues
s = 1;
a(1) = d(1);
P{1} = v(:,1)*v(:,1)';

for k = 2:length(d)
  if (abs(d(k) - d(k-1)) > qit.tol)
    % new eigenvalue, new projector
    s = s+1;
    a(s) = d(k);
    P{s} = v(:,k)*v(:,k)';
  else
    % extend current P
    P{s} = P{s} + v(:,k)*v(:,k)';
  end
end
