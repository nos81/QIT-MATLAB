function [A, spectrum] = liouvillian_fp(L, rtol)
% LIOUVILLIAN_FP  Fixed point states of a Liouvillian.
%  [A, spectrum] = liouvillian_fp(L, tol)
%
%  Finds the intersection of the kernel of the Liouvillian L
%  with the set of valid state operators, thus giving the set of
%  fixed point states for the quantum channel represented by the
%  master equation
%
%   \dot{\rho} = inv_vec(L * vec(\rho)).
%
%  Let size(L) == [D, D] (and d = sqrt(D) be the dimension of the Hilbert space).
%
%  Returns the D * n matrix A, which contains as its columns a set
%  of n vectorized orthogonal Hermitian matrices (with respect to
%  the Hilbert-Schmidt inner product) which "span" the set of FP
%  states in the following sense:
%
%    vec(\rho) = A * c,  where c \in \R^n and c_1 = 1.
%
%  A(:,1) is the shortest vector in the Hermitian kernel of L that
%  has trace 1, the other columns of A are traceless and normalized.
%  Hence, A defines an (n-1)-dimensional hyperplane.
%
%  A valid state operator also has to fulfill \rho \ge 0. These
%  operators form a convex set in the Hermitian trace-1 hyperplane
%  defined by A. Currently this function does nothing to enforce
%  positivity, it is up to the user to choose the coefficients a_k
%  such that this condition is satisfied.
%
%  Singular values of L less than or equal to the tolerance tol are
%  treated as zero.

% If L has Lindblad form, if L(rho) = \lambda * rho,
% we also have L(rho') = conj(\lambda) * rho'
% Hence the non-real eigenvalues of L come in conjugate pairs.
%
% Especially if rho \in Ker L, also rho' \in Ker L.

% Ville Bergholm 2011-2012


if nargin < 2
    rtol = eps(class(L));
end

% Hilbert space dimension
d = sqrt(size(L, 2));

% Get the kernel of L, restricted to vec-mapped Hermitian matrices.
% columns of A: complex orthogonal vectors A_i spanning the kernel (with real coefficients)
[A, spectrum] = nullspace_hermitian(L, rtol);

% Extract the trace-1 core of A and orthonormalize the rest.
% We want A_i to be orthonormal wrt. the H-S inner product
% <X, Y> := trace(X' * Y) = vec(X)' * vec(Y)

% With this inner product, trace(A) = <eye(d), A>.
temp = vec(speye(d));
a = (temp' * A)';

% Construct the shortest linear combination of A_i which has tr = 1.
% This is simply (1/d) * eye(d) _iff_ it belongs to span(A_i).
core = A * (a / norm(a)^2);  % trace-1

% Remove the component parallel to core from all A_i.
temp = core / norm(core); % normalized
A = A -temp * (temp' * A);

% Re-orthonormalize the vectors, add the core
A = [core, orthonormalize(A, 1e-6)]; % FIXME tolerance

N = size(A, 2);
for k = 1:N
  temp = inv_vec(A(:, k));
  A(:, k) = vec(0.5 * (temp + temp')); % re-Hermitize to fix numerical errors
end

% TODO intersect with positive ops
%probe(A);


%[u,t] = schur(L);
%unnormality = norm(t, 'fro')^2 -norm(diag(t))^2
%[us, ts] = ordschur(u, t, 'udi')
%E = ordeig(t);
% TODO eigendecomposition, find orthogonal complement to span(v) = ker(v').
% these are the states which do not belong to eigenspaces and
% should show transient polynomial behavior
end


function probe(A)
% Probe the positivity of the FP set

N = size(A, 2);

% upper limit for amplitudes of orthogonal non-core state components (from purity):
% necessary but not sufficient
C = A(:, 1);
a_max = sqrt(1 - C'*C);

x = 1.1 * a_max;
a = linspace(-x, x, 57);
b = linspace(-x, x, 52);
neg = [];
for j=1:length(b)
    for k=1:length(a)
        r = C +a(k)*A(:, 2) +b(j)*A(:, 3);
        neg(j, k) = min(eig(inv_vec(r)));
    end
end
figure();
mesh(a, b, neg);
xlabel('A_1')
ylabel('A_2')
zlabel('min eigenvalue')
end


function plot_eigenvalues(s)
% Plots the eigenvalues of L on the complex plane.

% propagate the eigenvalues a bit to see their behavior
t = linspace(0, 1, 10);
%for k = 1:length(t)
%  temp = expm(t(k) * L);
%  s_expL(:, k) = sort(eig(temp));
%end
s_expL = exp(s * t); % FIXME assuming the normality/diagonalizability of L...
if isreal(s_expL)
  s_expL = complex(s_expL);
end

figure();
% plot the unit circle
phi = linspace(0, 2*pi, 51);
plot(exp(1i * phi), '--');
%axis([-2 2 -2 2]);
hold on
plot(s, 'rx')
plot(s_expL.', 'o-')
legend('Unit circle', '\sigma(L)', '\sigma(exp(Lt))')
title('Eigenvalue plot')
axis equal
grid on
end
