function [S] = squeeze(z, n)
% SQUEEZE  Bosonic squeezing operator.
%  S = squeeze(z, n)
%
%  Returns the n-dimensional approximation for the bosonic
%  squeezing operator S(z) in the number basis {|0>, |1>, ..., |n-1>}.

% Ville Bergholm 2010


if (~isscalar(z))
  error('z must be a scalar.')
end

a = boson_ladder(n);
S = expm(0.5*(z'*a^2 -z*(a')^2));
