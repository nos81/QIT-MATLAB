function [C, s] = nullspace_hermitian(A, rtol)
% NULLSPACE_HERMITIAN  Kernel of a superoperator matrix restricted to the Hermitian subspace.
%  [C, spectrum] = nullspace_hermitian(A, rtol)
%
%  Solves the intersection of the kernel of the superop A and the Hermitian subspace.
%  A maps d*d matrices to whatever.
%
%  Singular values <= rtol * max(s) are considered zero.
    
% Hermitian and antihermitian (orthogonal) real subspaces: V = H \oplus A, h \in H, a \in A
% If G = A'*A is either block-diagonal or block-antidiagonal wrt these, 
% x = h+a, Gx = \lambda x implies that Gh = \lambda h and Aa = \lambda a.
% Hence...
% Ax := [Q, x] is block-antidiagonal if Q is hermitian
    
% Ville Bergholm 2011-2012


if nargin < 2
    rtol = eps(class(A));
end

D = size(A, 2);
d = sqrt(D); % since A is a superop


% find A restricted to the Hermitian (real) subspace H
V = reshape(H_basis(d), [D, D]); % Hermitian basis
U = [real(V); imag(V)]; % same, except now in a real vector space

%AH = U' * to_real(A) * U;
AH = to_real(A) * U;

% solve the kernel
[Z, s] = nullspace(AH, rtol); % null space (kernel) of AH

%spectrum_delta = s -svd(full(A)) % TODO for valid Lindblad ops these are equal, why?

% columns of C: complex orthogonal vectors spanning the kernel (with real coefficients)
C = V * Z;
end



function R = to_real(C)
% Represents a complex linear map in a real vector space.
% Corresponding transformation for vectors: x_R = [real(x); imag(x)];

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
