function P = projector(s)
% PROJECTOR  Projection operator defined by the state.
%  P = projector(s)
%
%  Returns the projection operator P defined by the state s.

% Ville Bergholm 2009-2014


global qit;

if (abs(purity(s) - 1) > qit.tol)
  error('The state is not pure, and thus does not correspond to a projector.')
end

s = to_op(s);
P = lmap(s);
