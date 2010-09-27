function U = phase(theta, dim)
% PHASE  Phase shift gate.
%  U = phase(theta)
%
%  Returns the phase shift gate U = diag(exp(i*theta)).

% Ville Bergholm 2010


n = length(theta);

if (nargin < 2)
  dim = n;
elseif (prod(dim) ~= n)
  error('Dimension mismatch.');
end

U = lmap(spdiags(exp(i*theta).', 0, n, n), {dim, dim});
