function D = trace_dist(r, s)
% STATE/TRACE_DIST  Trace distance of two states.
%  D = trace_dist(r, s)
%
%  Trace distance between state operators r and s is defined as
%  $D(r, s) = 0.5*\trace(\sqrt{A^\dagger * A})$, where $A = r-s$.
%
%  Equivalently $D(r, s) = 0.5*\sum_k |\lambda_k|$, where $\lambda_k$
%  are the eigenvalues of A (since A is Hermitian).

% Ville Bergholm 2009
%! M.A. Nielsen, I.L. Chuang, "Quantum Computation and Quantum Information" (2000), chapter 9.2.1


r = to_op(r);
s = to_op(s);

A = r.data -s.data;
D = 0.5*sum(abs(eig(A)));
%D = 0.5*trace(sqrtm(A'*A));
