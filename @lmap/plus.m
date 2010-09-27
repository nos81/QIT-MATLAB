function s = plus(s, t)
% PLUS  Addition of lmaps.
%  q = plus(s, t)
%
%  Returns the sum of lmaps s and t.

% Ville Bergholm 2010


if (~is_compatible(s, t))
  error('The lmaps are not compatible.');
end

s.data = s.data + t.data;
