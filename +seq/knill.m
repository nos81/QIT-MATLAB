function [s] = knill(phi)
% KNILL  Sequence for robust pi pulses.
%  s = knill(phi)
%
%  The target rotation in the NMR notation is \pi_\phi followed by Z_{-\pi/3}.
%  In an experimental setting the Z rotation can often be absorbed by a
%  reference frame change that does not affect the measurement results.
%
%  The Knill pulse is quite robust against off-resonance errors, and somewhat
%  robust against pulse strenght errors.

%! See the reference RHC2010.
% Ville Bergholm 2015-2016


if nargin < 1
  phi = 0; % default is R_x(pi)
end

theta = pi;

s = seq.nmr([theta, pi/6+phi;
             theta, phi;
             theta, pi/2+phi;
             theta, phi;
             theta, pi/6+phi]);
