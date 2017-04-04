function [W, a, b] = wigner(s, res, limits)
% WIGNER  Wigner quasi-probability distribution.
%  W = wigner(s, alpha)
%  [W, a, b] = wigner(s [, res, limits])
%
%  Returns the Wigner quasi-probability distribution
%  W(Im \alpha, Re \alpha) corresponding to the harmonic
%  oscillator state s given in the number basis.
%
%  The integral of W is normalized to unity.
%
%  NOTE: The truncation of the number state space to a finite dimension
%  results in spurious circular ripples in the Wigner function outside
%  a given radius. To increase the accuracy, increase the state space dimension.

% Ville Bergholm 2010


if nargin == 2 && isscalar(res)
  % return W(alpha)
  a = real(res);
  b = imag(res);
else
  % return a grid of W values
  if (nargin < 3)
    limits = [-2, 2, -2, 2];
    if nargin < 2
      res = [40, 40];
    end
  end

  a = linspace(limits(1), limits(2), res(1));
  b = linspace(limits(3), limits(4), res(2));
end

n = prod(s.dim);

% parity operator
P = ones(n,1);
P(2:2:n) = -1;
P = spdiags((2/pi)*P, 0, n, n); % include Wigner normalization here for convenience

W = zeros(length(b), length(a));
for k=1:length(a)
  for j=1:length(b)
    alpha = a(k)+i*b(j);
    r = s.prop(ho.displace(-alpha, n));
    %temp = r.data;
    %W(j,k) = real(temp'*(P.*temp));
    W(j,k) = real(r.ev(P));
  end
end
