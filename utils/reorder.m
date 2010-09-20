function [T, dim] = reorder(T, dim, p)
% REORDER  Change the order of subsystems in a tensor.
%  [T, dim] = reorder(T, dim, perm);
%    reorder(T, dim, [2 4 3 1]); % reorder the subsystems of T to the order given
%    reorder(T, dim, [2 5]);     % swap subsystems 2 and 5
%
%  Reorders the subsystems of the tensor T and the corresponding dimensions
%  in the vector dim according to permutation vector perm. All the tensor indices
%  are similarly reordered.
%  If only two subsystems are listed, swaps them.
%
%  Note: For now, only handles order-1 and order-2 tensors.

% Ville Bergholm 2009-2010


% number of systems
n = length(dim);

% total dimension
dd = prod(dim);

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
d = fliplr(dim);
dim = dim(perm); % now reorder the dimensions vector
perm = fliplr(n+1 -perm);

T = full(T); % FIXME Matlab reshape is broken (sparse matrices cause problems with singleton dimensions)

if (isvector(T))
  % vector
  % vector into tensor, reordering of dimensions, back into a "flat" vector
  T = reshape(permute(reshape(T, d), perm), [dd 1]);
else
  % assume matrix
  % TODO could also handle tensors of higher order
  % matrix into tensor, reordering of dimensions, back into a "flat" matrix
  T = reshape(permute(reshape(T, [d d]), [perm n+perm]), [dd dd]);
end
