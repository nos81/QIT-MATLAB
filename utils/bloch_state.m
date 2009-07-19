function s = bloch_state(a, dim)
% BLOCH_STATE  State corresponding to a Bloch vector.
%  s = bloch_state(a)      % assume dim == [sqrt(length(a) + 1)]
%  s = bloch_state(a, dim) % give state dimensions explicitly
%
%  Returns the state corresponding to the Bloch vector a.
%  The vector is defined in terms of the tensor basis
%  corresponding to dim.
%
%  For proper states norm(a) <= sqrt(prod(s.dim)-1) (e.g. for qubits norm(a) <= 1).

% Ville Bergholm 2009


n = length(a);

if (nargin == 1)
  dim = sqrt(n+1); % n == dim^2-1
end

G = tensorbasis(dim);
d = prod(dim);

a = a/sqrt(d); % to match the usual Bloch vector normalization
rho = eye(d)/d;
for k=1:n
  rho = rho +a(k)*G{k};
end
s = state(rho, dim);
