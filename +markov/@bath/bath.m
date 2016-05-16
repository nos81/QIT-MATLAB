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
  b.int_end = b.cut_omega * 500;

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



function plot_cut(b, boltz)
% Plots stuff as a function of cut limit
% boltz is the Boltzmann factor exp(-\beta \hbar \omega).

    om = -log(boltz) / b.scale  % \omega * TU

    % try different cutoffs relative to om
    N = 50;
    cut = linspace(0.5, 20, N);
    gs = zeros(2, N);
    gsm = zeros(2, N);
    for k=1:length(cut)
        k
        b.set_cutoff('smooth', cut(k) * abs(om));
        gs(:,k) = b.compute_gs(om);
        gsm(:,k) = b.compute_gs(-om);
    end
    odd_S = gs(2,:) -gsm(2,:);
    even_S = gs(2,:) +gsm(2,:);
    q = sign(om) ./ cut;
    figure();
    t = b.plot_stuff(cut, om, q, gs, odd_S, even_S);
    xlabel('cutoff relative to omega')
    title(sprintf('%s, boltzmann: %g, omega: %g', t, boltz, om));
end


function t = plot_stuff(b, x, om, q, gs, odd_S, even_S)

    boltz = exp(-b.scale * om);
    % ratio of Lamb shift to dephasing rate
    ratio = odd_S ./ (gs(1,:) .* (boltz+1));

    % gamma, S, ratio
    [ax, h1, h2] = plotyy(x, gs, x, ratio);
    hold(ax(1),'on');
    hold(ax(2),'on');
    set(h1,'Marker','s')
    set(h2,'Marker','x')

    % symmetrized and antisymmetrized S
    plot(ax(1), x, odd_S, 'k-', x, even_S, 'm-')

    % analytical expressions for even and odd S funcs
    switch b.cut_type
      case 'sharp'
        odd_S = log(abs(q.^2./(q.^2-1)));
        even_S = log(abs((1+q)./(1-q))) -2./q;
      case 'smooth'
        odd_S = 2*log(abs(q))./(1+q.^2);
        even_S = -pi./q  ./(1+q.^2);
      case 'exp'
        odd_S = -real(expint(q).*exp(q) +expint(-q).*exp(-q));
        even_S = real(expint(q).*exp(q) -expint(-q).*exp(-q)) -2./q;
      otherwise
        error('zzz')
    end
    plot(ax(1), x, om .* odd_S, 'ko', x, om .* even_S, 'mo');

    axis(ax, 'tight')
    set(get(ax(1),'Ylabel'),'String','[1/TU]')
    set(get(ax(2),'Ylabel'),'String','Lamb shift ratio')
    legend('gamma', 'S', 'S(\omega)-S(-\omega)', 'S(\omega)+S(-\omega)',...
           'S(\omega)-S(-\omega) (fermion)', 'S(\omega)+S(-\omega) (boson)', 'Lamb shift ratio')
    t = sprintf('Bath correlation tensor Gamma: %s, %s, relative T: %g', b.type, b.stat, 1/b.scale);
    grid(ax(1), 'on')
end


function build_LUT(b, om)
  % Build gs_table.

  % TODO justify limits for S lookup
  if nargin < 2
      % Default sampling for the lookup table.
      %lim = b.cut_omega;
      lim = log(10) / 5 / b.scale;  % up to boltzmann factor == 10
      om = logspace(log10(1.1 * lim), log10(5 * lim), 20); % logarithmic sampling
      om = [linspace(0.05 * lim, 1 * lim, 20), om]; % sampling is denser near zero, where S changes more rapidly
      om = [-fliplr(om), 0, om];  % symmetric around zero
  end

  b.omega = om;
  b.gs_table = zeros(2, length(om));
  for k=1:length(om)
    k
    b.gs_table(:,k) = b.compute_gs(om(k));
  end

  if 1
      % plot the LUT data
      S = b.gs_table(2,:);
      odd_S  = S -fliplr(S);
      even_S = S +fliplr(S);
      q = om / b.cut_omega;
      figure();
      t = b.plot_stuff(om, om, q, b.gs_table, odd_S, even_S);
      xlabel('omega [1/TU]')
      title(sprintf('%s, cutoff: %s, %g', t, b.cut_type, b.cut_omega));
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
