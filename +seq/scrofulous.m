function [s] = scrofulous(theta, phi)
% SCROFULOUS  Sequence for fixing errors in pulse length.
%  s = scrofulous(theta [, phi])
%
%  Returns the SCROFULOUS control sequence s for fixing errors
%  in pulse duration (or amplitude).
%
%  The target rotation is \theta_\phi in the NMR notation.
%
%  SCROFULOUS: Short Composite ROtation For Undoing Length Over- and UnderShoot

%! Cummins et al., "Tackling systematic errors in quantum logic gates with composite rotations", PRA 67, 042308 (2003).
% Ville Bergholm 2006-2011


if (nargin < 2)
  phi = 0; % default is R_x
end

th1 = fsolve(@(t) (sin(t)/t -(2/pi)*cos(theta/2)), 0.1);
ph1 = acos(-pi*cos(th1)/(2*th1*sin(theta/2)));
ph2 = ph1 - acos(-pi/(2*th1));

s = seq.nmr([th1, ph1; pi, ph2; th1, ph1]);
