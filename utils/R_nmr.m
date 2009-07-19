function R = R_nmr(theta, phi)
% R_NMR  SU(2) rotation \theta_\phi (NMR notation).
%  R = R_nmr(theta, phi)
%
%  Returns the one-qubit rotation by angle theta about the unit
%  vector [cos(phi), sin(phi), 0], or \theta_\phi in the NMR notation.

% Ville Bergholm 2009


global qit;

R = expm(-i * theta/2 * (cos(phi)*qit.sx +sin(phi)*qit.sy));
