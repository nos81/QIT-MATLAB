function U = two(B, t, dim)
% TWO  Two-qudit operator.
%
%  U = two(B, t, dim)
%
%  Returns the operator U corresponding to the two-qudit operator B applied
%  to subsystems t == [t1, t2] (and identity applied to the remaining subsystems).
%
%  dim is either the dimension vector for U or an integer scalar denoting
%  the number of subsystems in an all-qubit system.

% James Whitfield 2010
% Ville Bergholm 2010


n = length(dim);

if (length(t) ~= 2)
  error('Exactly two target subsystems required.')
end

if (any(t < 1) || any(t > n) || t(1) == t(2))
  error('Bad target subsystem(s).')
end

temp = prod(dim(t));

if (size(B,1) ~= temp || size(B,2) ~= temp)
  error('Dimensions of the target subsystems are not compatible with the dimension of U.')
end

a = min(min(t));
b = max(max(t));

if (t(1) < t(2))
  p = [1 3 2];
else
  p = [2 3 1];
end

% dimensions for t1, t2 and the subsystems in between
ddd = [dim(t), prod(dim(a+1:b-1))];

temp = reorder(kron(B, speye(ddd(3))), ddd, p);

U = mkron(speye(prod(dim(1:a-1))), temp, speye(prod(dim(b+1:end))));
