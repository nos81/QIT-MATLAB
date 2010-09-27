function s = mtimes(a, b)
% MTIMES  Multiplication of the state by a scalar.
%  q = times(a, s)
%  q = times(s, a)
%
%  Returns the state s multiplied by the scalar a.

% Ville Bergholm 2010


% a must be a state, otherwise we wouldn't be in this function
if (isa(b, 'lmap'))
  s = mtimes@lmap(a, b); % HACK for u_propagate, returns an lmap
elseif (isscalar(a) && isnumeric(a))
  s = b;
  s.data = a * s.data;
elseif (isscalar(b) && isnumeric(b))
  s = a;
  s.data = b * s.data;
else
  error('A state can only be multiplied by a scalar.')  
end
