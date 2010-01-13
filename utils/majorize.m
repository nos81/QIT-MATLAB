function res = majorize(x, y)
% MAJORIZE  Majorization partial order of real vectors.
%  res = majorize(x, y)
%
%  Returns true iff vector x is majorized by vector y.

% Ville Bergholm 2010


global qit

if (~isvector(x) || ~isvector(y))
  error('Both inputs must be vectors.')
end
    
if (length(x) ~= length(y))
  error('The vectors must be of equal length.')
end

x = cumsum(sort(x, 'descend'));
y = cumsum(sort(y, 'descend'));

if (abs(x(end) -y(end)) < qit.tol)
  % exact majorization
  res = all(x <= y);
else
  % weak majorization could still be possible, but...
  disp('Note: Vectors have unequal sums.')
  res = false;
end
