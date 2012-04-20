function plot_state_trajectory(traj, ls, reset)
% PLOT_STATE_TRAJECTORY  Plot a state trajectory in the correlation representation.
%  plot_state_trajectory(traj [, linestyle, reset])
%
%  For a single-qubit system, plots the trajectory in the Bloch sphere.
%
%  For a two-qubit system, plots the reduced single-qubit states (in
%  Bloch spheres), as well as the interqubit correlations.
%
%  traj is a cell vector of generalized Bloch vectors.
%  It can be obtained e.g. by using one of the continuous-time
%  state propagation functions and feeding the results to
%  bloch_vector.
%
%  If reset is false, adds another trajectory to current plot
%  without erasing it.
%
%  Example 1: trajectory of s under the Hamiltonian H
%   out = propagate(s, H, t, @(s,H) bloch_vector(s))
%   plot_bloch_trajectory(out)
%
%  Example 2: just a single state s
%   plot_bloch_trajectory({bloch_vector(s)})

% Ville Bergholm  2006-2012

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
    plot_traj(A, [2, 3, 4], ls);
  
  case 16
    % two qubits (or a single ququat...)
    if reset
      subplot(2,3,1)
      plot_bloch_sphere();
      title('qubit A')
    
      subplot(2,3,2)
      plot_bloch_sphere();
      title('qubit B')

      subplot(2,3,3)
      plot_2q_corr('diagonal');
      
      subplot(2,3,4)
      plot_2q_corr('pos');
      
      subplot(2,3,5)
      plot_2q_corr('neg');
    end

    subplot(2,3,1)
    plot_traj(A, [2, 3, 4], ls);

    subplot(2,3,2)
    plot_traj(A, [5, 9, 13], ls);

    subplot(2,3,3)
    plot_traj(A, [6, 11, 16], ls);

    subplot(2,3,4)
    plot_traj(A, [8, 10, 15], ls);

    subplot(2,3,5)
    plot_traj(A, [14, 7, 12], ls);
    
  otherwise
    error('At the moment only plots one- and two-qubit trajectories.')
end
end


function plot_traj(A, ind, ls)
% Plots the trajectory formed by the correlations given in ind.

  % if we only have a single point, do not bother with these
  if size(A, 2) > 1
      plot3(A(ind(1),1), A(ind(2),1), A(ind(3),1), [ls(1), 'x']);
      plot3(A(ind(1),:), A(ind(2),:), A(ind(3),:), ls);
      plot3(A(ind(1),end), A(ind(2),end), A(ind(3),end), [ls(1), 'o']);
  else
      plot3(A(ind(1),end), A(ind(2),end), A(ind(3),end), ls);
  end
end
