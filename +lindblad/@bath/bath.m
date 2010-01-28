classdef bath < handle
% BATH  Handle class for heat baths.
%
%  Currently only one type of bath is supported, a bosonic
%  canonical ensemble at absolute temperature T, with a
%  single-term coupling to the system.
%
%  The bath spectral density is Ohmic with a cutoff.
%    J(omega) = hbar^2*omega*cut(omega)*heaviside(omega);
%
%  Two types of cutoffs are supported: exponential and sharp.
%    cut(omega) = exp(-omega/omega_c)
%    cut(omega) = heaviside(omega_c - omega)

% gamma(omega) == 2*pi/hbar^2 (J(omega)-J(-omega))(1+n(omega))

% Ville Bergholm 2009-2010


properties (SetAccess = protected)
  % basic parameters
  type   % Bath type. Currently only 'ohmic' is supported.
  omega0 % hbar*omega0 is the unit of energy for all Hamiltonians (omega0 in Hz)
  T      % Absolute temperature of the bath (in K).

  % shorthand
  scale  % dimensionless temperature scaling parameter hbar*omega0/(kB * T)

  % additional parameters with default values
  cut_type  % spectral density cutoff type
  cut_limit % spectral density cutoff angular frequency/omega0
    
  % spectral density
  % J(omega0 * x)/omega0 = \hbar^2 * j(x) * heaviside(x) * cut_func(x);
  j        % spectral density profile
  cut_func % cutoff function

  % lookup tables / function handles for g and s
  % gamma(omega0 * x)/omega0 = g_func(x) * cut_func(x)
  % gamma(omega0 * dh(k))/omega0 = gs_table(1,k)
  %     S(omega0 * dh(k))/omega0 = gs_table(2,k)
  g_func
  g0
  s0
  dH
  gs_table
end  


methods
function b = bath(type, omega0, T)
  % constructor
  %  b = bath(type, omega0, T)
  %
  %  Sets up a descriptor for a heat bath coupled to a quantum system.
    
    
  global qit;

  if (nargin ~= 3)
    error('Usage: bath(type, omega0, T)')
  end

  % basic bath parameters
  b.type   = type;
  b.omega0 = omega0;
  b.T      = T;

  % shorthand
  b.scale = qit.hbar * omega0 / (qit.kB * T);

  switch type
    case 'ohmic'
      % Ohmic bosonic bath, canonical ensemble, single-term coupling

      b.g_func = @(x) 2*pi*x.*(1 +1./(exp(b.scale*x)-1));
      b.g0 = 2*pi/b.scale; % limit of g at x == 0
      b.j = @(x) x;

    otherwise
      error('Unknown bath type.')
  end

  % defaults, can be changed later
  set_cutoff(b, 'sharp', 20);
end


function build_LUT(b)
  % Build a lookup table for the S integral.

  error('unused')

  % TODO justify limits for S lookup
  if (limit > 10)
    temp = logspace(log10(10.2), log10(limit), 50);
    temp = [linspace(0.1, 10, 100), temp]; % sampling is denser near zero, where S changes more rapidly
  else
    temp = linspace(0.1, limit, 10);
  end
  b.dH = [-fliplr(temp), 0, temp]

  b.s_table = [];
  for k=1:length(b.dH)
    b.s_table(k) = S_func(b, b.dH(k));
  end

  b.dH      = [-inf, b.dH, inf];
  b.s_table = [0, b.s_table, 0];

  plot(b.dH, b.s_table, 'k-x');
end


function S = S_func(b, dH)
% Compute S(dH*omega0)/omega0 = P\int_0^\infty dv J(omega0*v)/(hbar^2*omega0) (dH*coth(scale/2 * v) +v)/(dH^2 -v^2)

  ep = 1e-5; % epsilon for Cauchy principal value
  if (abs(dH) <= 1e-8)
    S = b.s0;
  else
    % Cauchy principal value, integrand has simple poles at \nu = \pm dH.
    f = @(nu) b.j(nu) .* b.cut_func(nu) .* (dH*coth(b.scale * nu/2) +nu) ./ (dH^2 -nu.^2);
    S = quad(f, ep, abs(dH)-ep)...
        +quad(f, abs(dH)+ep, 100*b.cut_limit); % FIXME upper limit should be inf, this is arbitrary
  end
end


function set_cutoff(b, type, lim)
  b.cut_type = type;
  b.cut_limit = lim; % == omega_c/omega0

  % update cutoff function (at least Octave uses early binding, so when
  % parameters change we need to redefine it)
  switch b.cut_type
    case 'sharp'
      b.cut_func = @(x) (abs(x) <= b.cut_limit); % Heaviside theta cutoff
    case 'exp'
      b.cut_func = @(x) exp(-abs(x)/b.cut_limit); % exponential cutoff
    otherwise
      error('Unknown cutoff type "%s"', b.cut_type)
  end

  switch b.type
    case 'ohmic'
      b.s0 = -b.cut_limit; % limit of S at dH == 0
  end
  
  % clear lookup tables, since changing the cutoff requires recalc of S
  % start with a single point precomputed
  b.dH = [-inf, 0, inf];
  b.gs_table = [0, b.g0, 0;
		0, b.s0, 0];
end

end

end
