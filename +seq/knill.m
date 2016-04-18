function [s] = knill(theta, phi)
% KNILL  Robust pulses.
%  s = knill(theta [, phi])
%
%  The target rotation is \theta_\phi in the NMR notation.

%! TODO reference
% Ville Bergholm 2015-2016


if nargin < 2
  phi = 0; % default is R_x
end

s = seq.nmr([theta, pi/6+phi;
             theta, phi;
             theta, pi/2+phi;
             theta, phi;
             theta, pi/6+phi]);
