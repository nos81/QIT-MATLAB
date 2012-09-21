function s = u_propagate(s, P)
% U_PROPAGATE  Propagate the state a finite step in time using a propagator.
%  q = propagate(s, P)
%
%  Propagates the state s using the propagator P, which is either a
%  unitary Hilbert space propagator or a Liouville space propagator.
%  Returns the resulting state.

% Ville Bergholm 2009-2012

% TODO FIXME P as lmap
if isa(P, 'lmap')
  if is_ket(s)
    % state vector
    s = state(P * s);
  else
    % state operator
    s = state(P * s * P');
  end
  
  return
end

% assume P is a matrix

if length(P) == size(s.data, 1)
  % unitary
  if is_ket(s)
    % state vector
    s.data = P * s.data;
  else
    % state operator
    s.data = P * s.data * P';
  end
else
  % Liouvillian propagator
  s = to_op(s);
  s.data = inv_vec(P * vec(s.data));
end

