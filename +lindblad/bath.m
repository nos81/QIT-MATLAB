function [g, s] = bath(desc, dH)
% LINDBLAD/BATH  Bath response function. FIXME name
%
%  [g, s] = bath(bath_desc, dH)  % returns Gamma(omega0 * dH)/omega0, see below
%
%  Computes the bath response function Gamma for a single-term coupling to a
%  bosonic bath, a canonical ensemble at absolute temperature T. The bath has
%  an Ohmic spectral density with an exponential cutoff at omega_c:
%
%    J(omega) = hbar^2*omega*exp(-omega/omega_c)*heaviside(omega);
%
%  hbar*omega0 is the unit of energy for system Hamiltonians (omega0 in Hz).
%  It is stored in the bath descriptor along with other bath parameters.
%
%  When called with an initialized bath descriptor, it returns the corresponding
%  bath response function Gamma, evaluated at omega0 * dH:
%
%    Gamma(omega0 * dH)/omega0 == 0.5*g +i*s

% gamma(omega) == 2*pi/hbar^2 (J(omega)-J(-omega))(1+n(omega))
% == 2*pi*omega*exp(-abs(omega)/omega_c)(1+n(omega))

% Ville Bergholm 2009


if (nargin ~= 2)
  error('Need bath descriptor and dH.')
end

% assume parameters are set and lookup table computed
s = desc.N * interp1(desc.dH, desc.S, dH, 'linear', 0);
  
if (abs(dH) <= 1e-8)
  g = desc.N * 2*pi/desc.scale; % limit at dH == 0
else
  g = desc.N * 2*pi*dH*exp(-abs(dH)/desc.dH_c)*(1 +1/(exp(desc.scale*dH)-1));
end
