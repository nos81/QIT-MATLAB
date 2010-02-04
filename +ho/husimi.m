function [H, a, b] = husimi(s, varargin)
% HUSIMI  Husimi probability distribution.
%  H = husimi(s, alpha, z)
%  [H, a, b] = husimi(s [, res, limits, z])
%
%  Returns the Husimi probability distribution
%  H(Im \alpha, Re \alpha) corresponding to the harmonic
%  oscillator state s given in the number basis.
%
%  z is the optional squeezing parameter for the reference state:
%   |\alpha, z> = D(\alpha) S(z) |0>
%
%   H(s, \alpha, z) =  1/\pi <\alpha, z| \rho_s |\alpha, z>
%
%  The integral of H is normalized to unity.

% Ville Bergholm 2010


z = 0;
if (nargin == 3 && isscalar(varargin{1}) && isscalar(varargin{2}))
  % return H(alpha)
  alpha = varargin{1};
  a = real(alpha);
  b = imag(alpha);
  z = varargin{2};
else
  % return a grid of H values
  res = [40, 40];
  limits = [-2, 2, -2, 2];
  
  for k=1:(nargin-1)
    temp = varargin{k};
    
    switch (length(temp))
      case 1
        z = temp;
      case 2
        res = temp;
      case 4
        limits = temp;
      otherwise
        error('Unknown parameter type.')
    end
  end

  a = linspace(limits(1), limits(2), res(1));
  b = linspace(limits(3), limits(4), res(2));
end

n = prod(s.dim);
% normalization included here for convenience
ref = u_propagate(state(0, n)/sqrt(pi), ho.squeeze(z, n));

H = zeros(length(b), length(a));
for k=1:length(a)
  for j=1:length(b)
    alpha = a(k)+i*b(j);
    temp = u_propagate(ref, ho.displace(alpha, n));
    H(j,k) = fidelity(s, temp)^2;
  end
end
