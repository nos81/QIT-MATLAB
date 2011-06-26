function [s] = bb1(theta, phi, location)
% BB1  Sequence for fixing off-resonance (\sigma_z bias) errors.
%  s = bb1(theta [, phi])
%
%  Returns the Broadband number 1 control sequence s for fixing
%  errors in pulse lenght (or amplitude).
%
%  The target rotation is \theta_\phi in the NMR notation.

%! Cummins et al., "Tackling systematic errors in quantum logic gates with composite rotations", PRA 67, 042308 (2003).
% Ville Bergholm 2009

if (nargin < 3)
  location = 0.5; % default: symmetric
  if (nargin < 2)
    phi = 0; % default is R_x
  end
end

ph1 = acos(-theta/(4*pi));
W1  = seq.nmr([pi, ph1; 2*pi, 3*ph1; pi, ph1]);

%s = [W1; seq.nmr([theta, phi])];
s = [seq.nmr([location*theta, phi]); W1; seq.nmr([(1-location)*theta, phi])];
