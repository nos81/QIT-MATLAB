function ret = orthonormalize(A, atol)
% Orthonormalizes the basis in the matrix A (consisting of basis
% column vectors), possibly discarding linearly dependent ones.
%    
% Adaptation of the standard MATLAB function orth with variable tolerance.
    
% Ville Bergholm 2012


[U, S] = svd(A, 0);
[m, n] = size(A);
if m > 1
    s = diag(S);
elseif m == 1
    s = S(1);
else
    s = 0;
end
%tol = max(m,n) * max(s) * eps(class(A));
% We could use relative tolerance so scaling A doesn't change things,
% except it should. Almost-null vectors should be
% discarded, but if we have only those...

r = sum(s > atol);
ret = U(:, 1:r);
