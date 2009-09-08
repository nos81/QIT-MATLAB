function [desc, H, D] = bath_setup(type, varargin)
% LINDBLAD/BATH_SETUP  Create a bath descriptor.
%
%   bath_desc = bath_setup('ohmic', omega0, T, N [, omega_c = 10*omega0, max_dH = 10])
%
%   [bath_desc, H, D] = bath_setup('ohmic_fit',  Delta, T, T1, T2 [, omega_c = 10*omega0])
%
%  Sets up a descriptor for a bath coupled to a quantum system.
%
%  hbar*omega0 is the unit of energy for system Hamiltonians (omega0 in Hz), 
%  T (K) is the absolute temperature of the bath,
%  N is the dimensionless system-bath coupling factor squared,
%  omega_c (Hz) the cutoff angular frequency for the spectral density, and
%  max_dH the difference between the largest and smallest eigenvalue of the system Hamiltonian.
%
%  Alternatively, for a qubit system with decoherence times T1 and T2,
%  hbar*Delta is the system energy splitting. This form of the
%  function also returns the corresponding system Hamiltonian H and
%  the system-bath coupling operator D.
%
%  Currently only one type of bath is supported, a bosonic
%  canonical ensemble at absolute temperature T, with a
%  single-term coupling to the system. The bath spectral density is
%  Ohmic with an exponential cutoff at omega_c:
%
%    J(omega) = hbar^2*omega*exp(-omega/omega_c)*heaviside(omega);

% gamma(omega) == 2*pi/hbar^2 (J(omega)-J(-omega))(1+n(omega))
% == 2*pi*omega*exp(-abs(omega)/omega_c)(1+n(omega))

% Ville Bergholm 2009


global qit;

desc.type = type;
npar = length(varargin);

switch type
  case 'ohmic'
    % Ohmic bosonic bath, canonical ensemble, single-term coupling
    if (npar < 3)
      error('Required params: omega0, T, N.')
    end

    % define bath parameters
    desc.omega0 = varargin{1};
    desc.T      = varargin{2};
    desc.N      = varargin{3};
    
    desc.scale = qit.hbar * desc.omega0 / (qit.kB * desc.T); % physical temperature scaling parameter, dimensionless
    
    % optional params
    if (npar < 5)
      if (npar < 4)
        desc.omega_c = 10*desc.omega0;
      else
        desc.omega_c = varargin{4};
      end
    
      max_dH = 10;
    else
      max_dH = varargin{5};
    end
    desc.dH_c = desc.omega_c / desc.omega0;

    s0 = -desc.N * desc.dH_c; % limit of s at dH == 0
    J = @(x) x.*exp(-x/desc.dH_c);
    desc = build_LUT(desc, J, s0, max_dH);

    
  case 'ohmic_fit'
    % Fitting an ohmic bath to a given set of decoherence times

    if (npar < 4)
      error('Required params: Delta, T, T1, T2.')
    end

    % define bath parameters
    Delta = varargin{1}; % qubit energy splitting
    desc.omega0 = Delta/2;
    desc.T      = varargin{2};
    T1 = varargin{3};
    T2 = varargin{4};

    iTd = 1/T2 -0.5/T1; % inverse pure dephasing time
    if (iTd < 0)
      error('unphysical decoherence times')
    end
    
    desc.scale = qit.hbar * desc.omega0 / (qit.kB * desc.T); % physical temperature scaling parameter, dimensionless
    
    % optional params
    if (npar < 5)
      desc.omega_c = 10*desc.omega0;
    else
      desc.omega_c = varargin{5};
    end
    desc.dH_c = desc.omega_c / desc.omega0;

    % match bath couplings to T1, T2
    alpha = atan2(1, sqrt(T1*iTd*coth(desc.scale)*desc.scale*exp(-abs(Delta)/desc.omega_c)));
    desc.N = iTd*desc.scale/(desc.omega0*4*pi*cos(alpha)^2);

    % qubit Hamiltonian
    H = -qit.sz;
    max_dH = 10*desc.dH_c;

    % noise coupling
    D = cos(alpha)*qit.sz +sin(alpha)*qit.sx;

    % S lookup
    s0 = -desc.N * desc.dH_c; % limit of s at dH == 0
    J = @(x) x.*exp(-x/desc.dH_c);
    desc = build_LUT(desc, J, s0, max_dH);
    
    % decoherence times in scaled time units
    T1 = 1/(2*pi*desc.N*2* coth(desc.scale)* exp(-abs(Delta)/desc.omega_c) *sin(alpha)^2)
    T_dephase = desc.scale/(desc.N *4*pi*cos(alpha)^2);
    T2 = 1/(0.5/T1 +1/T_dephase)

  otherwise
    error('Unknown bath type.')
end

end



function desc = build_LUT(desc, J, s0, max_dH)
% Build a lookup table for the s integral.
% S(dH*omega0)/omega0 = P\int_0^\infty dv J(omega0*v)/(hbar^2*omega0) (dH*coth(scale/2 * v) +v)/(dH^2 -v^2)
%
% Function J(v) is actually J(omega0*v)/(hbar^2*omega0)

ep = 1e-5; % epsilon for Cauchy principal value

% TODO justify limits for S lookup
limit = max_dH * 1.05;
if (limit > 10)
  temp = logspace(log10(10.2), log10(limit), 50);
  desc.dH = [-fliplr(temp), linspace(-10, 10, 100), temp]; % sampling is denser near zero, where S changes more rapidly
else
  desc.dH = linspace(-limit, limit, 100);
end

desc.S = 0;
for k=1:length(desc.dH)
  dH = desc.dH(k);

  if (abs(dH) <= 1e-8)
    desc.S(k) = s0;
  else
    % Cauchy principal value, integrand has simple poles at +-dH.
    f = @(x) J(x).*(dH*coth(desc.scale*x/2) +x)./(dH^2 -x.^2);
    desc.S(k) = quad(f, ep, abs(dH)-ep)...
      +quad(f, abs(dH)+ep, 10*max_dH); % upper limit should be inf, this is arbitrary
  end
end

%plot(desc.dH, desc.S, 'k-x');
end
