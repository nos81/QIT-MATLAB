function s = mrdivide(s, a)
% MRDIVIDE  Division of the state by a scalar from the right.
%  q = rdivide(s, a)
%
%  Returns the state s divided by the scalar a.

% Ville Bergholm 2010


if (isscalar(a) && isnumeric(a))
  s.data = s.data / a;
else
  error('A state can only be divided by a scalar.')  
end
