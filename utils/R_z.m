function R = R_z(t)
% R_Z  SU(2) z-rotation.
%  U = R_z(theta)
%
%  Returns the one-qubit rotation about the z axis by the angle theta,
%  e^(-i \sigma_z theta/2).

% Ville Bergholm 2006-2009


R = [exp(-i*t/2), 0 ; 0, exp(i*t/2)];
