function s = u_propagate(s, U)
% U_PROPAGATE  Propagate the state using a unitary.
%  q = propagate(s, U)
%
%  Propagates the state s using the unitary propagator U,
%  returns the resulting state.

% Ville Bergholm 2009-2010


if (isa(U, 'lmap'))
  if (size(s.data, 2) == 1)
    % state vector
    s = state(U * s);
  else
    % state operator
    s = state(U * s * U');
  end
  return
end

% assume U is a matrix
if (size(s.data, 2) == 1)
  % state vector
  s.data = U*s.data;
else
  % state operator
  s.data = U*s.data*U';
end
