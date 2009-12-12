classdef bath
% BATH  Class for heat baths.
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

% Ville Bergholm 2009


properties (SetAccess = protected)
  % basic parameters
  type   % Bath type. Currently only 'ohmic' is supported.
  omega0 % hbar*omega0 is the unit of energy for all Hamiltonians (omega0 in Hz)
  T      % Absolute temperature of the bath (in K).

  % shorthand
  scale  % physical temperature scaling parameter, dimensionless
end

properties
  % additional parameters with default values
  cut_type  % spectral density cutoff type
  cut_limit % spectral density cutoff angular frequency (Hz)
  N         % dimensionless system-bath coupling factor squared
end

properties (SetAccess = protected)
  % spectral density
  % J(omega0 * x)/omega0 = N * \hbar^2 * j(x) * heaviside(x) * cut_func(x);
  j        % spectral density profile
  cut_func % cutoff function

  % lookup tables / function handles for g and s
  % gamma(omega0 * x)/omega0 = N * g_func(x) * cut_func(x)
  % S(omega0 * x)/omega0     = N * s_lut(x)
  g_func
  g0
  dH
  s_lut
  s0
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

  % defaults, can be changed later
  b.N = 0; % bath coupling strength
  b.cut_type  = 'sharp';
  b.cut_limit = 20;

  switch type
    case 'ohmic'
      % Ohmic bosonic bath, canonical ensemble, single-term coupling

      b.g_func = @(x) 2*pi*x*(1 +1/(exp(b.scale*x)-1));
      b.g0 = 2*pi/b.scale; % limit of g at x == 0
      b.s0 = -b.cut_limit; % limit of s at x == 0
      b.j = @(x) x;

    otherwise
      error('Unknown bath type.')
  end

  b = build_LUT(b, 20); % FIXME max_dH
end


function b = build_LUT(b, limit)
  % Build a lookup table for the S integral.
  % S(dH*omega0)/omega0 = P\int_0^\infty dv J(omega0*v)/(hbar^2*omega0) (dH*coth(scale/2 * v) +v)/(dH^2 -v^2)

  ep = 1e-5; % epsilon for Cauchy principal value

  % TODO justify limits for S lookup
  if (limit > 10)
    temp = logspace(log10(10.2), log10(limit), 50);
    b.dH = [-fliplr(temp), linspace(-10, 10, 100), temp]; % sampling is denser near zero, where S changes more rapidly
  else
    b.dH = linspace(-limit, limit, 100);
  end

  b.s_lut = [];
  for k=1:length(b.dH)
    dH = b.dH(k);

    if (abs(dH) <= 1e-8)
      b.s_lut(k) = b.s0;
    else
      % Cauchy principal value, integrand has simple poles at +-dH.
      f = @(x) b.j(x) .* b.cut_func(x) .* (dH*coth(b.scale*x/2) +x) ./ (dH^2 -x.^2);
      b.s_lut(k) = quad(f, ep, abs(dH)-ep)...
        +quad(f, abs(dH)+ep, 10*b.cut_limit); % FIXME upper limit should be inf, this is arbitrary
    end
  end

  plot(b.dH, b.s_lut, 'k-x');
end



function b = set.cut_limit(b, val)
  b.cut_limit = val; % == omega_c/omega0
  b = cut_update(b);
end


function b = set.cut_type(b, val)
  b.cut_type = val;
  b = cut_update(b);
end


function b = cut_update(b)
% update cutoff function (at least Octave uses early binding, so when parameters change we need to redefine it)

  switch b.cut_type
    case 'sharp'
      b.cut_func = @(x) (abs(x) <= b.cut_limit); % Heaviside theta cutoff
    case 'exp'
      b.cut_func = @(x) exp(-abs(x)/b.cut_limit); % exponential cutoff
    otherwise
      error('Unknown cutoff type "%s"', b.cut_type)
  end
  
  % clear lookup tables, since changing the cutoff requires recalc of S
  b.dH = [];
  b.s_lut  = [];
end

end

end
