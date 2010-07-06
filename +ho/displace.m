function [D] = displace(alpha, n)
% DISPLACE  Bosonic displacement operator.
%  D = displace(alpha, n)
%
%  Returns the n-dimensional approximation for the bosonic
%  displacement operator D(alpha) in the number basis {|0>, |1>, ..., |n-1>}.

% Ville Bergholm 2010


if (~isscalar(alpha))
  error('Alpha must be a scalar.')
end

a = boson_ladder(n);
D = expm(alpha*a' -alpha'*a);
