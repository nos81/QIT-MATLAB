function s = nmr(a)
% NMR  Convert NMR-style rotations into a one-qubit control sequence.
%  s = NMR([theta1, phi1; ...])
%
%  Returns a one-qubit control sequence corresponding to NMR rotations
%  of the form \theta_\phi.

%  [a, theta] ^= R_a(theta) = expm(-i*a*sigma*theta/2) = expm(-i*H*t) => H = a*sigma/2, t = theta

% Ville Bergholm 2006-2011


global qit;

i_theta = 1;
i_phi = 2;

% find theta angles that are negative, convert them to corresponding positive rotation
rows = find(a(:,i_theta) < 0);
a(rows,i_theta) = -a(rows,i_theta);
a(rows,i_phi) = a(rows,i_phi)+pi;

s.A = 0;
s.B = {1i*0.5*qit.sx, 1i*0.5*qit.sy};
s.tau = a(:, i_theta);
s.control = [cos(a(:, i_phi)), sin(a(:, i_phi))];
