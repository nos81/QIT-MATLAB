function [A, spectrum] = superop_fp(L, rtol)
% SUPEROP_FP  Fixed point states of a Liouvillian superoperator.
%  [A, spectrum] = superop_fp(L, tol)
%
%  Finds the intersection of the kernel of the Liouvillian L
%  with the set of valid state operators, thus giving the set of
%  fixed point states for the quantum channel represented by the
%  master equation
%
%   \dot{\rho} = inv_vec(L * vec(\rho)).
%
%  L maps the space of vectorized complex d \times d matrices to itself.
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
%  positivity, it is up to the user to choose the coefficients c_k
%  such that this condition is satisfied.
%
%  Singular values of L less than or equal to the tolerance tol are
%  treated as zero.

% If L has Lindblad form, if L(rho) = \lambda * rho,
% we also have L(rho') = conj(\lambda) * rho'
% Hence the non-real eigenvalues of L come in conjugate pairs.
%
% Especially if rho \in Ker L, also rho' \in Ker L.

% Ville Bergholm 2011-2016


if nargin < 2
    rtol = eps(class(L));
end

% Liouville space dimension
D = size(L, 2);
% Hilbert space dimension
d = sqrt(D);


%% Find the kernel of L, restricted to vec-mapped hermitian matrices.
% Since the hermitian matrices do not form a complex subspace of
% the complex liouville space but rather a real subspace, we must
% first map all the objects to a corresponding real vector space.

% columns of V: hermitian basis within the domain of L, vec-mapped, orthonormal
V = reshape(H_basis(d), [D, D]);  % V is unitary

% Find the restriction of L to the real hermitian subspace within its domain.
LH = vec_to_real(L*V);  % == map_to_real(L) * vec_to_real(V);

% solve the kernel of LH
% columns of KH: real orthonormal vectors spanning the kernel (with real coefficients)
[KH, spectrum] = nullspace(LH, rtol);

%spectrum_delta = svd(KH) -svd(L) % TODO for valid Lindblad ops these are equal, why?

% Map kernel back to a complex vector space. Columns of K are orthonormal.
K = V * KH;


%% Find the intersection of the kernel span_R(K) and the trace(X) == 1 hyperplane.
% Extract the shortest trace-1 vector of span_R(K) and
% orthonormalize the rest, using _real_ coefficients to preserve hermiticity!

% <X, Y> := trace(X' * Y) = vec(X)' * vec(Y)
% trace(X) = <I, X> = vec(I)' * vec(X) = 1. Normalizing, we obtain
% u' * vec(x) = 1/sqrt(d), where the unit vector u = vec(I)/sqrt(d).

% Project u, the normal vector of the tr=1 plane onto span(K), never mind normalization.
c = vec(speye(d));
% Coordinate vector w = K' * c  is always real since the coords are traces of hermitian matrices.
% cK = K * inv(K'*K) * K' * c
% Columns of K are orthonormal, hence K'*K = I.
cK = K * K' * c;
% scale so that trace(cK) = 1
cK = cK / (c' * cK);

% Remove the component parallel to cK from the kernel.
temp = cK / norm(cK); % unit vector
K = K -temp * (temp' * K);

% Re-orthonormalize the vectors, add the core
A = [cK, orthonormalize(K, 1e-6)]; % FIXME tolerance

test_H(A);
return

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


function test_H(A)
% tests the hermitianness of each column of A
N = size(A, 2);
for k = 1:N
  temp = inv_vec(A(:, k));
  if norm(temp-temp') > 1e-6
      disp(sprintf('col %d is non-hermitian', k))
  end
end
end


function R = vec_to_real(x)
% Represents a complex vector in a real vector space.
    R = [real(x); imag(x)];
end

function R = map_to_real(C)
% Represents a complex linear map in a real vector space.
    R = [real(C), -imag(C); imag(C), real(C)];
end


function U = H_basis(d)
% Basis for the Hermitian subspace within the vector space of complex d*d matrices.
% The antihermitian basis is obtained by multiplying the hermitian one with 1i.

D = d^2;

x = 1/sqrt(2);
U = zeros(d, d, D);
n = 1;

% loop over the lower triangle (incl. diagonal)
for b = 1:d
  for a = b:d
    if a == b
      % diagonals: real => hermitian
      U(a, b, n) = 1;
      n = n+1;
    else
      % offdiagonals
      % real, symmetric or imag, antisymmetric => hermitian
      U(a, b, n) = x;
      U(b, a, n) = x;
      n = n+1;
      U(a, b, n) = 1i*x;
      U(b, a, n) = -1i*x;
      n = n+1;
    end
  end
end
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
