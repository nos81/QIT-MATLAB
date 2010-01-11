function s = minus(s, t)
% STATE/MINUS  Subtraction of states.
%  q = minus(s, t)
%
%  Returns the difference of s and t.
%  The states are not normalized at any point.

% Ville Bergholm 2010


if any(s.dim ~= t.dim)
  error('The dimensions of the states do not match.')
end

if any(size(s.data) ~= size(t.data))
  error('Addition of pure and mixed states is not defined.')
end

s.data = s.data - t.data;
