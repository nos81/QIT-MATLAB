function R = R_y(t)
% R_Y  SU(2) y-rotation.
%  U = R_y(theta)
%
%  Returns the one-qubit rotation about the y axis by the angle theta,
%  e^(-i \sigma_y theta/2).

% Ville Bergholm 2006-2009


sy = [0 -i; i 0];
R = expm(-i*t/2*sy);
