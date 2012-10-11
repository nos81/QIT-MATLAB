function [s] = bb1(theta, phi, location)
% BB1  Sequence for correcting pulse length errors.
%  s = bb1(theta [, phi])
%
%  Returns the Broadband number 1 control sequence s for fixing
%  proportional errors in pulse lenght (or amplitude).
%
%  The target rotation is \theta_\phi in the NMR notation.

%! Cummins et al., "Tackling systematic errors in quantum logic gates with composite rotations", PRA 67, 042308 (2003).
% Ville Bergholm 2009-2011

if (nargin < 3)
  location = 0.5; % default: symmetric
  if (nargin < 2)
    phi = 0; % default is R_x
  end
end

ph1 = acos(-theta/(4*pi));
s  = seq.nmr([location*theta, phi; pi, ph1; 2*pi, 3*ph1; pi, ph1; (1-location)*theta, phi]);
