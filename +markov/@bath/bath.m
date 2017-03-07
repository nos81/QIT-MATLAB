classdef bath < handle
% BATH  Handle class for heat baths.
%
%  Currently supports bosonic and fermionic canonical ensembles at
%  absolute temperature T, with an ohmic spectral density.
%  The system-bath coupling is assumed to be of the form
%    H_int = A \otimes \sum_k \lambda_k (a_k +a_k')
%
%  Below, dimensional quantities are denoted using a backslash: \t = t * TU etc.
%  TU denotes the time unit used.
%
%  The bath spectral density is ohmic with a cutoff:
%    \J(\omega) = \omega * cut_func(\omega) * heaviside(\omega)
%
%  Three types of cutoffs are supported:
%    cut_func(x * \omega_c) =
%       exponential:  exp(-x)
%       smooth:       (1+x^2)^{-1}
%       sharp:        heaviside(1-x)
%
%  We represent the spectral correlation tensor \Gamma of the bath as follows.
%  Since we only have one coupling term, \Gamma is a scalar. It depends on three parameters:
%  the inverse temperature of the bath \s = \beta \hbar,
%  the spectral cutoff frequency of the bath \omega_c,
%  and the system frequency \omega. \Gamma has the following scaling property:
%
%    \Gamma_{\s,\omega_c}(\omega) = \Gamma_{\s/a,\omega_c*a}(\omega*a)/a.
%
%  Hence we may eliminate the dimensions by choosing a = TU:
%
%    \Gamma_{\s,\omega_c}(\omega) = 1/TU * Gamma_{s, omega_c}(omega).
%
%  Furthermore, we split Gamma to its hermitian and antihermitian parts (notice the prefactors):
%    Gamma(omega) = 0.5*gamma(omega) +1i*S(omega)

% Ville Bergholm 2009-2017


properties %(SetAccess = protected)
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

  pf        % Planck or Fermi function n(omega) for the bath, depending on bath statistics
  corr_int  % integral transform for the dimensionless bath correlation function

  % lookup tables / function handles for g and s
  % \gamma(omega(k)/TU) * TU = gs_table(1,k)
  %     \S(omega(k)/TU) * TU = gs_table(2,k)
  omega
  g_func
  s_func
  g0      % limit of g at x == 0
  s0      % limit of S at x == 0
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
    
  if nargin ~= 4
    error('Usage: bath(...)')
  end

  % basic bath parameters
  b.type   = type;
  b.stat   = stat;
  b.TU     = TU;
  b.T      = T;

  % dimensionless physical scale factor
  global qit;
  b.scale = qit.hbar / (qit.kB * T * TU);

  % bath spectral density
  switch b.type
    case 'ohmic'
    otherwise
      error('Unknown bath type.')
  end

  % defaults, can be changed later
  b.set_cutoff('exp', 1);
end


function set_cutoff(b, type, lim)
% Define the cutoff function for J.
% We assume that cut_func(0) == 1.
% If type or lim is [], keep the old value.

  if ~isempty(type)
    b.cut_type = type;
  end
  if ~isempty(lim)
      b.cut_omega = lim;  % == \omega_c * TU
  end

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
      b.pf = @(x) 1./(exp(b.scale * x) - 1);
      b.corr_int = @(s, nu) nu .* b.cut_func(nu) .* (exp(-1i*nu.*s) +2*cos(nu.*s).*b.pf(nu));
      b.g_func = @(x) 2*pi * x .* b.cut_func(abs(x)) .* (1 +b.pf(x));
      % s_func has simple poles at \nu = \pm x.
      %b.s_func = @(x,nu) nu .* b.cut_func(nu) .* ((1+b.pf(nu))./(x-nu) +b.pf(nu)./(x+nu));
      b.s_func = @(x,nu) nu .* b.cut_func(nu) .* (x * coth(b.scale * nu/2) +nu) ./ (x^2 -nu.^2);
      b.g0 = 2*pi / b.scale;
      b.s0 = -quad(b.cut_func, 0, b.int_end);

    case 'fermion'
      b.pf = @(x) 1./(exp(b.scale * x) + 1);
      b.corr_int = @(s, nu) nu .* b.cut_func(nu) .* (exp(-1i*nu.*s) +2i*sin(nu.*s).*b.pf(nu));
      b.g_func = @(x) 2*pi * abs(x) .* b.cut_func(abs(x)) .* (1 -b.pf(x));
      % s_func has simple poles at \nu = \pm x.
      %b.s_func = @(x,nu) nu .* b.cut_func(nu) .* ((1-b.pf(nu))./(x-nu) +b.pf(nu)./(x+nu));
      b.s_func = @(x,nu) nu .* b.cut_func(nu) .* (x +nu .* tanh(b.scale * nu/2)) ./ (x^2 -nu.^2);
      b.g0 = 0;
      b.s0 = -quad(@(x) b.cut_func(x) .* tanh(x*b.scale/2), 0, b.int_end);

    otherwise
      error('Unknown statistics.')
  end
  
  % clear lookup tables, initialize with the (known) limits at infinity and zero
  b.omega   = [-Inf, 0, Inf];
  b.gs_table = [[0; 0], [b.g0; b.s0], [0; 0]];
end


function t = desc(b)
% Bath description string.
    t = sprintf('%s, %s, beta hbar omega_c: %g, cutoff: %s, %g',...
                b.type, b.stat, b.scale * b.cut_omega, b.cut_type, b.cut_omega);
end


function res = plot_correlation(b)
% Plots the bath correlation function \C_{\s,\omega_c}(\t) = <B(\t) B(0)> / \hbar^2.
% It scales as
%   \C_{\s,\omega_c}(\t) = \C_{\s/a,\omega_c*a}(\t/a) / a^2.
% Choosing a = TU, we obtain
%   \C_{\s,\omega_c}(\t) = 1/TU^2 * C_{s, omega_c}(t).

    tol_nu = 1e-7;  % approaching the singularity at nu=0 this closely
    c = b.cut_omega;
    figure();

    % plot the functions to be transformed
    subplot(1,3,1)
    nu = linspace(tol_nu, 5*c, 500);
    plot(nu, nu .* b.cut_func(nu) .* b.pf(nu), 'r')
    hold on;
    if strcmp(b.stat, 'boson')
        plot(nu, nu .* b.cut_func(nu) .* (1+b.pf(nu)), 'b')
    else
        plot(nu, nu .* b.cut_func(nu) .* (1-b.pf(nu)), 'g')
    end
    grid on;
    legend('n', '1 \pm n')
    xlabel('$\nu$ [1/TU]')
    title('Integrand without exponentials')

    % plot the full integrand
    subplot(1,3,2)
    res = [];
    t = linspace(0, 4/c, 5);
    nu = linspace(tol_nu, 5*c, 100);
    for k=1:length(t)
        k
        fun = @(x) b.corr_int(t(k), x);
        res(:,k) = fun(nu);
    end
    plot(nu, real(res), '-', nu, imag(res), '--')
    grid on;
    xlabel('$\nu$ [1/TU]')
    title('Integrand')

    % plot the correlation function
    subplot(1,3,3)
    res = [];
    t = linspace(0, 10/c, 100);  % real part of C(t) is even, imaginary part odd
    for k=1:length(t)
        k
        %fun = @(nu) nu .* b.cut_func(nu) .* (exp(-1i*nu.*t(k)).*(1+b.pf(nu))+exp(1i*nu.*t(k)).*b.pf(nu));
        fun = @(x) b.corr_int(t(k), x);
        res(k,1) = quad(fun, tol_nu, b.int_end);
        %res(k,1) = quad(fun, tol_nu, 100);
    end
    plot(t, real(res), 'k-', t, imag(res), 'k--') %, t, abs(res), '-.')
    a = axis();
    hold on;
    grid on;
    xlabel('t [TU]')
    ylabel('[1/TU^2]')
    title(['Bath correlation function C(t), ', b.desc()])

    % plot analytic high- and low-temp limits
    x = c * t;
    switch b.cut_type
      case 'sharp'
        temp = ((1 +1i*x).*exp(-1i*x) -1)./x.^2;
        hb = 2/b.scale * sin(x)./x;
      case 'smooth'
        temp = -1i*sinh(x).*sinint(1i*x) +cosh(x).*(1i*pi/2+(expint(x)+expint(-x))/2)...
               -1i * pi/2 * exp(-abs(x));  % imag part
        hb = pi/b.scale * exp(-abs(x));
      case 'exp'
        temp = 1./(1+(1i*x)).^2;
        hb = 2/b.scale ./ (1+x.^2);
      otherwise
        error('Unknown cutoff type "%s"', b.cut_type)
    end
    temp = temp * c^2;
    hb   = hb * c;
    switch b.stat
      case 'boson'
        plot(t, real(temp), 'b.')  % real part, cold
        plot(t, hb, 'r.')  % real part, hot
        plot(t, imag(temp), 'ko')  % imag part, every T
        legend('re C', 'im C', 're C (cold, analytic)', 're C (hot, analytic)', 'im C (analytic)')
      case 'fermion'
        plot(t, real(temp), 'k.')  % real part, every T
        plot(t, imag(temp), 'bo')  % imag part, cold
        plot(t, 0, 'ro')  % imag part, hot
        legend('re C', 'im C', 're C (analytic)', 'im C (cold, analytic)', 'im C (hot, analytic)')
      otherwise
        error('Unknown statistics.')
    end
    axis(a)
    return

    temp = t(t < 1/c);
    plot(temp, res(1)*(1-(temp*c).^2), 'k-')
    temp = t(t >= 1/c & t < b.scale);
    plot(temp, res(1)*(temp*c).^(-1), 'k:')
    temp = t(t > b.scale);
    plot(temp, res(1)*exp(-temp./b.scale), 'k--')
end


function plot_LUT(b)
% Plots the LUT contents plus something extra as a function of omega.

    if length(b.omega) <= 3
        b.build_LUT();
    end

    om = b.omega;
    % computed values
    S = b.gs_table(2,:);
    odd_S  = S -fliplr(S);
    even_S = S +fliplr(S);

    figure();
    b.plot_stuff(om, om, b.gs_table, odd_S, even_S);
    hold on;
    % omega cutoff value as vertical lines
    a = axis();
    c = b.cut_omega;
    plot([c,c], a(3:4), 'k-')
    plot([-c,-c], a(3:4), 'k-')
    xlabel('omega [1/TU]')
    title(b.desc())
end


function plot_cutoff(b, boltz)
% Plots stuff as a function of the cutoff frequency.
% boltz is the Boltzmann factor exp(-\beta \hbar \omega) that fixes
% omega when bath temperature is fixed.

    orig_scale = b.scale;
    orig_cut   = b.cut_omega;

    om = -log(boltz) / b.scale  % \omega * TU
    % using the scaling property, scale with |om|
    temp = abs(om);
    om = om/temp;  % sign of om
    b.scale = orig_scale * temp;

    % try different cutoffs relative to om
    N = 50;
    cut = linspace(0.5, 5, N);  % cut_omega/|omega|
    gs = zeros(2, N);
    gsm = zeros(2, N);
    for k=1:length(cut)
        k
        b.set_cutoff([], cut(k));

        gs(:,k) = b.compute_gs(om);
        gsm(:,k) = b.compute_gs(-om);
    end
    gs  = gs  * temp;
    gsm = gsm * temp;
    odd_S = gs(2,:) -gsm(2,:);
    even_S = gs(2,:) +gsm(2,:);

    % restore original bath parameters
    b.scale = orig_scale;
    b.set_cutoff([], orig_cut);

    figure();
    b.plot_stuff(cut, om*ones(size(cut)), gs, odd_S, even_S);
    xlabel('\omega_c/\omega')
    t = b.desc();
    title(sprintf('%s, boltzmann: %g, omega: %g', t, boltz, om*temp));
end


function plot_stuff(b, x, om, gs, odd_S, even_S)
% Plot stuff as a function of the vector x. The other inputs are given as a function of x.

    boltz = exp(-b.scale * om);
    % ratio of Lamb shift to dephasing rate, Ising spin chain
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
    q = om / b.cut_omega;
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
    grid(ax(1), 'on')
end


function build_LUT(b, om)
  % Build gs_table.

  % TODO justify limits for S lookup
  if nargin < 2
      % Default sampling for the lookup table, up to 5*cut_omega
      lim = b.cut_omega;
      %lim = log(10) / 5 / b.scale;  % up to boltzmann factor == 10
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
