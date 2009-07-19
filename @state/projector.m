function P = projector(s)
% PROJECTOR  Projection operator defined by state.
%  P = projector(s)
%
%  Returns the projection operator P defined by the state s.

% Ville Bergholm 2009


if (size(s.data, 2) == 1)
  % state vector
  P = s.data * s.data';
else
  % state operator
  error('not a pure state')
  %P = s.data;
end
