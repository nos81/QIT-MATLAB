function plot_bloch_trajectory(traj, c)
% PLOT_BLOCH_TRAJECTORY  Plot a trajectory in a Bloch sphere.
%  plot_bloch_trajectory(trajectory, linestyle)
%
%  Plots the first three non-identity components of the trajectory
%  of a quantum system, obtained using one of the continuous-time
%  state propagation functions and feeding the results to
%  bloch_vector.
%
%  A Bloch sphere can be plotted separately using plot_bloch_sphere.
%
%  Example:
%   [out] = propagate(s, H, t, @bloch_vector)
%   plot_bloch_trajectory(out, linestyle)

% Ville Bergholm  2006-2011


if nargin < 2
    c = 'b-';
end

a = cell2mat(traj);

% discard identity component
plot3(a(2,:), a(3,:), a(4,:), c);
plot3(a(2,end), a(3,end), a(4,end), 'k.');
