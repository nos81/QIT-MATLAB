function res = majorize(x, y)
% MAJORIZE  Majorization partial order of real vectors.
%  res = majorize(x, y)
%
%  Returns true iff vector x is majorized by vector y,
%  i.e. res = x \preceq y.


% Ville Bergholm 2010-2012


global qit

if (~isvector(x) || ~isvector(y) || ~isreal(x) || ~isreal(y))
  error('Inputs must be real vectors.')
end

if (length(x) ~= length(y))
  error('The vectors must be of equal length.')
end

x = cumsum(sort(x, 'descend'));
y = cumsum(sort(y, 'descend'));

if (abs(x(end) -y(end)) <= qit.tol)
  % exact majorization
  res = all(x - y <= qit.tol);
else
  % weak majorization could still be possible, but...
  disp('Note: Vectors have unequal sums.')
  res = false;
end
