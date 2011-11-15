function plot_state_trajectory(traj, ls, reset)
% PLOT_STATE_TRAJECTORY  Plot a state trajectory in the correlation representation.
%  plot_state_trajectory(traj [, linestyle, reset])
%
%  For a single-qubit system, plots the trajectory in the Bloch sphere.
%
%  For a two-qubit system, plots the reduced single-qubit states (in
%  Bloch spheres), as well as the diagonal correlations.
%
%  traj is a cell vector of generalized Bloch vectors.
%  It can be obtained e.g. by using one of the continuous-time
%  state propagation functions and feeding the results to
%  bloch_vector.
%
%  If reset is false, add another trajectory to current plot
%  without erasing it.
%
%  Example:
%   out = propagate(s, H, t, @(s,H) bloch_vector(s))
%   plot_bloch_trajectory(out, linestyle)

% Ville Bergholm  2006-2011

if nargin < 3
    reset = true;

    if nargin < 2
      ls = 'b-';
    end
end

A = cell2mat(traj);

switch size(A, 1)
  case 4
    % single qubit
    if reset
      plot_bloch_sphere();
    end
    plot3(A(2,:), A(3,:), A(4,:), ls);
    plot3(A(2,end), A(3,end), A(4,end), 'k.');
  
  case 16
    % two qubits (or a single ququat...)
    if reset
      subplot(1,3,1)
      plot_bloch_sphere();
      title('qubit A')
    
      subplot(1,3,2)
      plot_bloch_sphere();
      title('qubit B')

      subplot(1,3,3)
      plot_2q_corr();
      title('diagonal correlations')
    end

    subplot(1,3,1)
    plot3(A(2,:), A(3,:), A(4,:), ls);
    plot3(A(2,end), A(3,end), A(4,end), 'k.');

    subplot(1,3,2)
    plot3(A(5,:), A(9,:), A(13,:), ls);
    plot3(A(5,end), A(9,end), A(13,end), 'k.');

    subplot(1,3,3)
    plot3(A(6,:), A(11,:), A(16,:), ls);
    plot3(A(6,end), A(11,end), A(16,end), 'k.');
    
  otherwise
    error('At the moment only plots one- and two-qubit trajectories.')
end
