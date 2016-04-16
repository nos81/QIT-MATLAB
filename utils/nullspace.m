function [Z, s] = nullspace(A, rtol)
% NULLSPACE  Solves the kernel of a matrix.
%  [Z, s] = nullspace(A, rtol)
%
%  Given a matrix A and a relative tolerance rtol, returns a basis
%  for the kernel (null space) of A in the columns of Z, and
%  optionally the singular values of A in s.
%
%  Adaptation of the standard MATLAB function null, with variable tolerance.

% Ville Bergholm 2011-2016


if nargin < 2
    rtol = eps(class(A));
end

[m, n] = size(A);
[~, S, V] = svd(A, 0);
if m > 1
    s = diag(S);
elseif m == 1
    s = S(1);
else
    s = 0;
end
%tol = max(m,n) * max(s) * eps(class(A));
% use relative tolerance so scaling A doesn't change things
% cond = max(s) / min(s);  % condition number of A TODO

tol = rtol * max(s);
r = sum(s > tol);
Z = V(:, r+1:n);
