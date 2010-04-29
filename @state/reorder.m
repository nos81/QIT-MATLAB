function s = reorder(s, p)
% REORDER  Change the order of subsystems in a state.
%  x = reorder(s, perm);
%  x = reorder(s, [2 4 3 1]); % reorder the subsystems of s to the order given
%  x = reorder(s, [2 5]);     % swap subsystems 2 and 5
%
%  Reorders the subsystems of state s according to permutation vector perm.
%  If only two subsystems are listed, swaps them.

% Ville Bergholm 2009


if (nargin < 2)
  error('Need the permutation.')
end

% number of systems
n = length(s.dim);
% total dimension
dd = prod(s.dim);

perm = 1:n;
if (length(p) == 2)
  % swap
  perm(p(1)) = p(2);
  perm(p(2)) = p(1);
else
  if (length(setxor(perm, p)) ~= 0)
    error('Invalid permutation.')
  end
  perm = p;
end

% big-endian ordering is more natural for users, but matlab funcs
% prefer little-endian, so we reverse it
d = fliplr(s.dim);
s.dim = s.dim(perm); % now reorder the dimensions vector
perm = fliplr(n+1 -perm);

if (size(s.data, 2) == 1)
  % state vec
  % flat vector into tensor, reordering of dimensions, back into a flat vector
  s.data = reshape(permute(reshape(s.data, d), perm), [dd 1]);
else
  % state op
  % flat matrix into tensor, reordering of dimensions, back into a flat matrix
  s.data = reshape(permute(reshape(s.data, [d d]), [perm n+perm]), [dd dd]);
end
