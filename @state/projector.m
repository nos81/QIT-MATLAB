function P = projector(s)
% PROJECTOR  Projection operator defined by the state.
%  P = projector(s)
%
%  Returns the projection operator P defined by the state s.

% Ville Bergholm 2009-2010


global qit;

if (abs(purity(s) - 1) > qit.tol)
  error('The state is not pure, and thus does not correspond to a projector.')
end

if is_ket(s)
  % state vector
  P = s.data * s.data';
else
  % state operator
  P = s.data;
end
