function s = ptranspose(s, sys)
% PTRANSPOSE  Partial transpose.
%  x = ptranspose(s, sys);
%  x = ptranspose(s, [1 3]);  % transpose subsystems 1 and 3 of state s
%
%  Returns the partial transpose of the state s 
%  wrt. the subsystems listed in the vector sys.

% TODO not sure if this should return a state object, since a ptransposed state may no longer be a valid quantum state.

% Ville Bergholm 2008-2010


if (nargin < 2)
  error('Need the partitioning.')
end

% number of systems
dim = dims(s);
n = length(dim);

% total dimension
dd = prod(dim);

% which systems to transpose, into binary vector
tran = zeros(1,n);
tran(sys) = 1;

% big-endian ordering is more natural
d = fliplr(dim);
tran = fliplr(tran);

% swap the transposed dimensions
for k=1:n
  if (tran(k))
    % transposed
    perm(k) = n+k;
    perm(n+k) = k;
  else
    % not transposed
    perm(k) = k;
    perm(n+k) = n+k;
  end
end

s = to_op(s);

% flat matrix into tensor, partial transpose, back into a flat matrix
s.data = reshape(permute(reshape(s.data, [d d]), perm), [dd dd]);
