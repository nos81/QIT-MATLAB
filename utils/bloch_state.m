function s = bloch_state(a, dim)
% BLOCH_STATE  State corresponding to a generalized Bloch vector.
%  s = bloch_state(A)      % assume dim == sqrt(size(A))
%  s = bloch_state(A, dim) % give state dimensions explicitly
%
%  Returns the state s corresponding to the generalized Bloch vector A.
%
%  The vector is defined in terms of the standard Hermitian tensor basis B
%  corresponding to the dimension vector dim.
%
%    \rho_s == \sum_{ijk...} A_{ijk...} B_{ijk...} / \sqrt(D),
%
%  where D = prod(dim). For valid states norm(A) <= sqrt(D).

% Ville Bergholm 2009-2011

s = size(a);
if (isvector(a))
  s = length(a); % HACK, we don't want singleton dims.
end
n = prod(s);

if (nargin == 1)
  dim = sqrt(s); % s == dim.^2
end

G = tensorbasis(dim);
d = prod(dim);

a = a/sqrt(d); % to match the usual Bloch vector normalization
rho = zeros(d);
for k=1:n
  rho = rho +a(k)*G{k};
end
s = state(rho, dim);
