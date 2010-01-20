function [W, a, b] = wigner(s, res, limits)
% WIGNER  Wigner quasi-probability distribution.
%  W = wigner(s, alpha)
%  [W, a, b] = wigner(s [, res, limits])
%
%  Returns the Wigner quasi-probability distribution
%  W(Im \alpha, Re \alpha) corresponding to the harmonic
%  oscillator state s.
%
%  Example: pcolor(wigner(state(0, 20)))

% Ville Bergholm 2010


if (nargin == 2 && isscalar(res))
  % return W(alpha)
  a = real(res);
  b = imag(res);
else
  % return a grid of W values
  if (nargin < 3)
    limits = [-2, 2, -2, 2];
    if (nargin < 2)
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
    r = u_propagate(s, ho.displace(-alpha, n));
    %temp = r.data;
    %W(j,k) = real(temp'*(P.*temp));
    W(j,k) = real(ev(r, P));
  end
end

%pcolor(a, b, W);
%shading interp;
%shading flat;
%xlabel('Re(\alpha)')
%ylabel('Im(\alpha)')
