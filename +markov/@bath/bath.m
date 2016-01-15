classdef bath < handle
% BATH  Handle class for heat baths.
%
%  Currently supports bosonic and fermionic canonical ensembles at
%  absolute temperature T, with an ohmic spectral density.
%
%  The system-bath coupling is assumed to be of the form
%    H_int = A \otimes \sum_k \lambda_k (a_k +a_k')
%
%  The bath spectral density is ohmic with a cutoff:
%    J(\omega) = \omega * cut_func(\omega) * heaviside(\omega)
%
%  Three types of cutoffs are supported:
%    cut_func(x * \omega_c) =
%       exponential:  exp(-x)
%       smooth:       (1+x^2)^{-1}
%       sharp:        heaviside(1-x)

% \gamma(omega) == 2*pi (J(omega)-J(-omega))(1+n(omega))

% Ville Bergholm 2009-2016


properties (SetAccess = protected)
  % basic parameters
  type   % Bath type. Currently only 'ohmic' is supported.
  stat   % Bath statistics ('boson' or 'fermion').
  TU     % Time unit (in s). All Hamiltonians have been made dimensionless by multiplying with TU/\hbar.
  T      % Absolute temperature of the bath (in K).

  % shorthand
  scale  % dimensionless temperature scaling parameter \hbar/(kB * T * TU)

  % additional parameters with default values
  cut_type  % spectral density cutoff type
  cut_omega % spectral density cutoff angular frequency (in 1/TU)
    
  % spectral density
  % J(omega/TU) * TU = omega * heaviside(omega) * cut_func(omega);
  cut_func % cutoff function

  % lookup tables / function handles for g and s
  % \gamma(omega/TU) * TU = g_func(omega)
  %      S(omega/TU) * TU = s_func(omega)
  % \gamma(omega(k)/TU) * TU = gs_table(1,k)
  %      S(omega(k)/TU) * TU = gs_table(2,k)
  g_func
  s_func
  g0      % limit of g at x == 0
  s0      % limit of S at x == 0
  omega
  gs_table
  
  % HACKs
  int_end  % integration interval endpoint for improper integrals
end  

properties (SetAccess = public)
  int_pv = 1;
end


methods (Static)
function y = interpolate(ee, tt, x)
  % interp1 does way too many checks
  y = tt(:,1) + ((x - ee(1))/(ee(2) - ee(1)))*(tt(:,2) - tt(:,1));
end
end


methods
function b = bath(type, stat, TU, T)
  % constructor
  %  b = bath(type, TU, T)
  %
  %  Sets up a descriptor for a heat bath coupled to a quantum system.
    
    
  global qit;

  if nargin ~= 4
    error('Usage: bath(...)')
  end

  % basic bath parameters
  b.type   = type;
  b.stat   = stat;
  b.TU     = TU;
  b.T      = T;

  % dimensionless physical scale factor
  b.scale = qit.hbar / (qit.kB * T * TU);

  % bath spectral density
  switch b.type
    case 'ohmic'
    otherwise
      error('Unknown bath type.')
  end

  % defaults, can be changed later
  b.set_cutoff('exp', 10);
end


function set_cutoff(b, type, lim)
% Define the cutoff function for J
% We assume that cut_func(0) == 1.

  b.cut_type = type;
  b.cut_omega = lim;  % == \omega_c * TU

  % TODO quad cannot handle improper integrals, switch to int() or quadgk()
  b.int_end = b.cut_omega * 100;

  % update cutoff function (Matlab uses early binding, so when parameters change we need to redefine it)
  switch b.cut_type
    case 'sharp'
      b.cut_func = @(x) heaviside(b.cut_omega -x); % Heaviside theta cutoff
      b.int_end = b.cut_omega; % no need to integrate further
    case 'smooth'
      b.cut_func = @(x) 1./(1 +(x/b.cut_omega).^2);  % rational cutoff
    case 'exp'
      b.cut_func = @(x) exp(-x / b.cut_omega);     % exponential cutoff
    otherwise
      error('Unknown cutoff type "%s"', b.cut_type)
  end
  b.setup();
end


function setup(b)
% Initialize the g and s functions, and the LUT.
% Must be called after parameters change.

  % bath statistics
  switch b.stat
    case 'boson'
      b.g_func = @(x) 2*pi * x .* b.cut_func(abs(x)) .* (1 +1./(exp(b.scale * x)-1));
      % s_func has simple poles at \nu = \pm x.
      %b.s_func = @(x,nu) nu .* b.cut_func(nu) .* ((1+temp(nu))./(x-nu) +temp(nu)./(x+nu));
      b.s_func = @(x,nu) nu .* b.cut_func(nu) .* (x * coth(b.scale * nu/2) +nu) ./ (x^2 -nu.^2);
      b.g0 = 2*pi / b.scale;
      b.s0 = -quad(b.cut_func, 0, b.int_end);
      
    case 'fermion'
      temp = @(x) 1./(exp(b.scale * x) + 1);
      b.g_func = @(x) 2*pi * abs(x) .* b.cut_func(abs(x)) .* (1 -temp(x));
      % s_func has simple poles at \nu = \pm x.
      %b.s_func = @(x,nu) nu .* b.cut_func(nu) .* ((1-temp(nu))./(x-nu) +temp(nu)./(x+nu));
      b.s_func = @(x,nu) nu .* b.cut_func(nu) .* (x +nu .* tanh(b.scale * nu/2)) ./ (x^2 -nu.^2);
      b.g0 = 0;
      b.s0 = -quad(@(x) b.cut_func(x) .* tanh(x*b.scale/2), 0, b.int_end);
      
    otherwise
      error('Unknown statistics.')
  end
  
  % clear lookup tables
  % start with a single point precomputed
  b.omega = [-Inf, 0, Inf];
  b.gs_table = [0, b.g0, 0;
		0, b.s0, 0];
end


function build_LUT(b, omegas)
  % Build gs_table.

  % TODO justify limits for S lookup
  if nargin < 2
    omegas = logspace(log10(1.01 * b.cut_omega), log10(6 * b.cut_omega), 20);
    omegas = [linspace(0.1, 0.99 * b.cut_omega, 20), omegas]; % sampling is denser near zero, where S changes more rapidly
  end
  b.omega = [-fliplr(omegas), 0, omegas];
  om = b.omega;  % shorthand

  b.gs_table = zeros(2, length(om));
  for k=1:length(om)
    k
    b.gs_table(:,k) = b.compute_gs(om(k));
  end

  if true
      plot(om, b.gs_table, '-o');
      hold on
      temp = b.gs_table(2,:);
      delta_S = temp-fliplr(temp);
      plot(om, delta_S, 'k-')
      % analytical \Delta S
      q = om / b.cut_omega;
      switch b.cut_type
        case 'sharp'
          temp = log(abs(q.^2./(q.^2-1)));
        case 'smooth'
          temp = 2*log(abs(q))./(1+q.^2);
        case 'exp'
          temp = -real(expint(q).*exp(q) +expint(-q).*exp(-q));
        otherwise
          error('zzz')
      end
      plot(om, om .* temp, 'ko');
      xlabel('omega [1/TU]')
      ylabel('[1/TU]')
      legend('gamma', 'S', 'S(\omega)-S(-\omega)')
      title(sprintf('Bath correlation tensor Gamma: %s, %s, cutoff: %s, %g, relative T: %g', b.type, b.stat, b.cut_type, b.cut_omega, 1/b.scale));
      grid on
  end

  % add the limits at infinity
  b.omega   = [-Inf, b.omega, Inf];
  b.gs_table = [[0; 0], b.gs_table, [0; 0]];
end


function ret = compute_gs(b, omega)
% Compute and return [\gamma(omega/TU) * TU, S(omega/TU) * TU]

    ep = 1e-7; % epsilon for Cauchy principal value
    tol_omega0 = 1e-8;

    ret = [0; 0];

    if abs(omega) <= tol_omega0
      ret = [b.g0; b.s0];  % limit at omega == 0
    else
      ret(1) = b.g_func(omega);

      % Cauchy principal value. Integrand has simple poles at \nu = \pm omega,
      % only the positive one hits the integration region.
      fun = @(nu) b.s_func(omega, nu);
      if b.int_pv
          ret(2) = quad(fun, tol_omega0, abs(omega)-ep)...
                   +quad(fun, abs(omega)+ep, b.int_end);
          % TODO justify tol_omega0 here
      else
          % Complex path integration method
          d = 5;
          p1 = quad(fun, 0, 1i*d);
          p2 = quad(fun, 1i*d, b.int_end+1i*d);
          ret(2) = real(p1)+real(p2);
      end
    end
end


function [g, s] = corr(b, x)
% CORR  Bath spectral correlation tensor.
%
%  [g, s] = bath.corr(omega)
%
%  Returns the bath spectral correlation tensor evaluated at omega,
%  split into the real and imaginary parts:
%
%    \Gamma(omega / TU) * TU == 0.5*g +1i*s

  tol_omega = 1e-8;
  max_ip_omega = 0.1; % maximum interpolation distance TODO justify


  % assume parameters are set and lookup table computed
  %s = interp1(b.dH, b.s_table, x, 'linear', 0);

  % use the lookup table
  % binary search for the interval in which omega falls
  first = 1;
  last = length(b.omega);
  while first+1 ~= last
    pivot = round((first + last)/2);
    if x < b.omega(pivot)
      last = pivot;
    else
        first = pivot;
    end
  end
  ee = b.omega([first, last]);
  tt = b.gs_table(:, [first, last]);
  % now x is in [ee(1), ee(2))

  % distances of x to the endpoints
  gap = ee(2) -ee(1);
  d1 = abs(x -ee(1));
  d2 = abs(x -ee(2));

  % do we need to compute a new point, or just interpolate?
  if d1 <= tol_omega  % close enough to either endpoint?
    temp = b.gs_table(:, first);

  elseif d2 <= tol_omega
    temp = b.gs_table(:, last);

  elseif gap <= max_ip_omega + tol_omega  % short enough gap to interpolate?
    temp = b.interpolate(ee, tt, x);

  else  % compute a new point, add it into the LUT
    if gap <= 2*max_ip_omega  % split the gap in half
      p = ee(1) +gap/2;
      if x < p
        idx = 2; % which ee p will replace
      else
          idx = 1;
      end
    elseif d1 <= max_ip_omega  % max interpolation distance from endpoint
      p = ee(1)+max_ip_omega;
      idx = 2;
    elseif d2 <= max_ip_omega
      p = ee(2)-max_ip_omega;
      idx = 1;
    else  % just use x as the new point
      p = x;
      idx = 1;
    end

    % compute new g, s values at p
    temp = b.compute_gs(p);

    % insert them into the LUT
    b.omega = [b.omega(1:first), p, b.omega(last:end)];
    b.gs_table = [b.gs_table(:, 1:first), temp, b.gs_table(:, last:end)];

    % now interpolate the required value
    ee(idx) = p;
    tt(:, idx) = temp;
    temp = b.interpolate(ee, tt, x);
  end

  g = temp(1);
  s = temp(2);
end

end
end
