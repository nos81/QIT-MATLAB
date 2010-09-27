function U = two(B, t, dim)
% TWO  Two-qudit operator.
%
%  U = two(B, t, dim)
%
%  Returns the operator U corresponding to the bipartite operator B applied
%  to subsystems t == [t1, t2] (and identity applied to the remaining subsystems).
%
%  dim is the dimension vector for U.

% James Whitfield 2010
% Ville Bergholm 2010


n = length(dim);

if (length(t) ~= 2)
  error('Exactly two target subsystems required.')
end

if (any(t < 1) || any(t > n) || t(1) == t(2))
  error('Bad target subsystem(s).')
end

dB = B.dim;

temp = dim(t);
if (~isequal(dB{2}, temp))
  error('Dimensions of the target subsystems are not compatible with the dimensions of B.')
end

a = min(min(t));
b = max(max(t));

% how tensor(B_12, I_3) should be reordered
if (t(1) < t(2))
  p = [1 3 2];
else
  p = [2 3 1];
end

% dimensions for the untouched subsystems
before    = prod(dim(1:a-1));
inbetween = prod(dim(a+1:b-1));
after     = prod(dim(b+1:end));

B = reorder(tensor(B, lmap(speye(inbetween))), {p, p});

U = tensor(lmap(speye(before)), B, lmap(speye(after)));

% restore dimensions
d1 = dim;
d1(t) = dB{1};
U = lmap(U, {d1, dim});
