function U = two(B, t, d_in)
% TWO  Two-qudit operator.
%
%  U = two(B, t, d_in)
%
%  Returns the operator U corresponding to the bipartite-to-bipartite operator B applied
%  to subsystems t == [t1, t2] (and identity applied to the remaining subsystems).
%
%  d_in is the input dimension vector for U.

% James Whitfield 2010
% Ville Bergholm 2010-2011


n = length(d_in);

if (length(t) ~= 2)
  error('Exactly two target subsystems required.')
end

if (any(t < 1) || any(t > n) || t(1) == t(2))
  error('Bad target subsystem(s).')
end

dB = B.dim;

if ~isequal(dB{2}, d_in(t))
  error('Dimensions of the target subsystems are not compatible with the input dimensions of B.')
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
before    = prod(d_in(1:a-1));
inbetween = prod(d_in(a+1:b-1));
after     = prod(d_in(b+1:end));

% tensor in the corresponding identities
B = reorder(tensor(B, lmap(speye(inbetween))), {p, p});
U = tensor(lmap(speye(before)), B, lmap(speye(after)));

% restore dimensions
d_out = d_in;
d_out(t) = dB{1};
U = lmap(U, {d_out, d_in});
