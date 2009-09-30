function b = bath(type, omega0, T)
% LINDBLAD/BATH  Create a bath object.
%
%   b = bath('ohmic', omega0, T)
%
%  Sets up a descriptor for a heat bath coupled to a quantum system.
%
%  hbar*omega0 is the unit of energy for system Hamiltonians (omega0 in Hz), 
%  T is the absolute temperature of the bath (in K).
%
%  Additional properties of the bath can be set afterwards using the set method.
%  These include
%
%  N: dimensionless system-bath coupling factor squared,
%  cut_type:
%  cut_limit:  (Hz) the cutoff angular frequency for the spectral density, and
%  max_dH the difference between the largest and smallest eigenvalue of the system Hamiltonian.
%
%  Currently only one type of bath is supported, a bosonic
%  canonical ensemble at absolute temperature T, with a
%  single-term coupling to the system. The bath spectral density is
%  Ohmic with an exponential cutoff at omega_c:
%
%    J(omega) = hbar^2*omega*exp(-omega/omega_c)*heaviside(omega);

% gamma(omega) == 2*pi/hbar^2 (J(omega)-J(-omega))(1+n(omega))
% == 2*pi*omega*exp(-abs(omega)/omega_c)(1+n(omega))

% TODO bath should really be a class... too bad they don't work in octave without hacks
% Ville Bergholm 2009


global qit;

if (nargin ~= 3)
  error('Usage: bath(type, omega0, T)')
end

% basic bath parameters
b.type   = type;
b.omega0 = omega0;
b.T      = T;

% shorthand
b.scale = qit.hbar * omega0 / (qit.kB * T); % physical temperature scaling parameter, dimensionless


% defaults, can be changed later
b.N = 0; % bath coupling strength
b   = lindblad.bath_set(b,  'cut_type', 'sharp',  'cut_limit', 20);
% FIXME changing the cutoff requires recalc of S

switch type
  case 'ohmic'
    % Ohmic bosonic bath, canonical ensemble, single-term coupling

    b.g_func = @(x) 2*pi*x*(1 +1/(exp(b.scale*x)-1));
    b.g0 = 2*pi/b.scale; % limit of g at x == 0

    s0 = -b.cut_limit; % limit of s at x == 0
    J = @(x) x .* b.cut_func(x); % J(omega0 * x)/omega0

  otherwise
    error('Unknown bath type.')
end

b = build_LUT(b, J, s0, 10); % FIXME max_dH
end



function b = build_LUT(b, J, s0, max_dH)
% Build a lookup table for the s integral.
% S(dH*omega0)/omega0 = P\int_0^\infty dv J(omega0*v)/(hbar^2*omega0) (dH*coth(scale/2 * v) +v)/(dH^2 -v^2)
%
% Function J(v) is actually J(omega0*v)/(hbar^2*omega0)

ep = 1e-5; % epsilon for Cauchy principal value

% TODO justify limits for S lookup
limit = max_dH;
if (limit > 10)
  temp = logspace(log10(10.2), log10(limit), 50);
  b.dH = [-fliplr(temp), linspace(-10, 10, 100), temp]; % sampling is denser near zero, where S changes more rapidly
else
  b.dH = linspace(-limit, limit, 100);
end

b.S = [];
for k=1:length(b.dH)
  dH = b.dH(k);

  if (abs(dH) <= 1e-8)
    b.S(k) = s0;
  else
    % Cauchy principal value, integrand has simple poles at +-dH.
    f = @(x) J(x) .* (dH*coth(b.scale*x/2) +x) ./ (dH^2 -x.^2);
    b.S(k) = quad(f, ep, abs(dH)-ep)...
      +quad(f, abs(dH)+ep, 5*b.cut_limit); % FIXME upper limit should be inf, this is arbitrary
  end
end

%plot(b.dH, b.S, 'k-x');
end
