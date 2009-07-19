function [s] = corpse(theta, phi)
% CORPSE  Sequence for fixing off-resonance (\sigma_z bias) errors.
%  s = corpse(theta [, phi])
%
%  Returns the CORPSE control sequence s for fixing off-resonance
%  errors, i.e. ones arising from a constant but unknown
%  \sigma_z bias in the Hamiltonian.
%
%  The target rotation is \theta_\phi in the NMR notation.
%
%  CORPSE: Compensation for Off-Resonance with a Pulse SEquence

%! Cummins et al., "Tackling systematic errors in quantum logic gates with composite rotations", PRA 67, 042308 (2003).
% Ville Bergholm 2009


if (nargin < 2)
  phi = 0; % default is R_x
end

n = [1 1 0]; % CORPSE
%n = [0 1 0]; % short CORPSE

temp = asin(sin(theta/2)/2);

th1 = 2*pi*n(1) +theta/2 -temp;
th2 = 2*pi*n(2) -2*temp;
th3 = 2*pi*n(3) +theta/2 -temp;

s = seq.nmr([th1, phi;
	 th2, phi+pi;
	 th3, phi]);
