function [g, s] = bath_corr(b, dH)
% LINDBLAD/BATH_CORR  Bath spectral correlation tensor.
%
%  [g, s] = bath_corr(bath, dH)  % returns Gamma(omega0 * dH)/omega0, see below
%
%  Returns the bath spectral correlation tensor Gamma evaluated at omega0 * dH:
%
%    Gamma(omega0 * dH)/omega0 == 0.5*g +i*s

% gamma(omega) == 2*pi/hbar^2 (J(omega)-J(-omega))(1+n(omega))

% Ville Bergholm 2009


if (nargin ~= 2)
  error('Need a bath object and dH.')
end

% assume parameters are set and lookup table computed
s = b.N * interp1(b.dH, b.S, dH, 'linear', 0);
  
if (abs(dH) <= 1e-8)
  g = b.N * b.g0; % limit at dH == 0
else
  g = b.N * b.g_func(dH) * b.cut_func(dH);
end
