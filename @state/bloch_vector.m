function a = bloch_vector(s, as_tensor)
% BLOCH_VECTOR  Generalized Bloch vector.
%  A = bloch_vector(s)
%
%  Returns the generalized Bloch vector A corresponding to the state s.
%
%  For an n-subsystem state the generalized Bloch vector is an order-n correlation
%  tensor defined in terms of the standard Hermitian tensor basis B
%  corresponding to s.dim:
%
%    A_{ijk...} == \sqrt(D) * \trace(\rho_s  B_{ijk...}),
%
%  where D = prod(s.dim). A is always real since \rho_s is Hermitian.
%  For valid states norm(A) <= sqrt(D) (e.g. for a qubit system norm(A) <= 2).

% Ville Bergholm 2009-2011


dim = dims(s);
G = tensorbasis(dim);
n = length(G);

a = zeros(n, 1);
for k=1:n
  a(k) = ev(s, G{k});
end
a = a * sqrt(prod(dim)); % to match the usual Bloch vector normalization

% should we return it as a vector, or as a tensor?
if nargin < 2 || ~as_tensor
  return
end

% into an array, one dimension per subsystem
if (length(dim) == 1)
  dim = [dim, 1]; % stupid reshape input syntax exceptions
end
a = reshape(a, dim.^2);
