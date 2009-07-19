function R = R_x(t)
% R_X  SU(2) x-rotation.
%  U = R_x(theta)
%
%  Returns the one-qubit rotation about the x axis by the angle theta,
%  e^(-i \sigma_x theta/2).

% Ville Bergholm 2006-2009


sx = [0 1; 1 0];
R = expm(-i*t/2*sx);
