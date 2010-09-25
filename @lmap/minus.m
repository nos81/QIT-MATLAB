function s = minus(s, t)
% MINUS  Subtraction of lmaps.
%  q = minus(s, t)
%
%  Returns the difference of s and t.
%  The lmaps are not normalized at any point.

% Ville Bergholm 2010


if (~is_compatible(s, t))
  error('The lmaps are not compatible.');
end

s.data = s.data - t.data;
