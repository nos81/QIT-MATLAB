function [h, ind] = plot_2q_corr(labels)
% PLOT_2Q_CORR  Plots the correlations simplex for two-qubit states.
%  [h, ind] = plot_2q_corr(group)
%
%  Plots the geometrical representation of the set of allowed
%  correlations in a two-qubit state. For each group of three
%  correlation variables the set is a tetrahedron.
%
%  The groups are 'diagonal', 'pos' and 'neg'.
%  For diagonal correlations the vertices correspond to the four Bell states.
%
%  Returns the graphics handle h and a vector of three linear
%  indices denoting the correlations to be plotted as the x, y and z coordinates.

% NOTE the strange logic in the ordering of the pos and neg
% correlations follows the logic of the Bell state labeling convention, kind of.

% Ville Bergholm  2011-2012


if nargin < 1
  labels = 'diagonal';
end

% vertices and faces
v = [1, 1, -1;  1, -1, 1;  -1, 1, 1; -1, -1, -1];
f = [1,2,4; 2,3,4; 3,1,4; 3,2,1];

hold off;
% tetrahedron
h = patch('Faces', f, 'Vertices', v, 'FaceVertexCData', 0, 'FaceColor', 'flat');
hold on;
alpha(0.2)
grid on
view(-35, 30)
axis equal


% mark vertices
plot3(0, 0, 0, 'r.'); % center
plot3(v(:,1), v(:,2), v(:,3), 'r.');

% label axes and vertices
switch labels
  case 'diagonal'
    title('diagonal correlations')
    xlabel('XX');
    ylabel('YY');
    zlabel('ZZ');
    text(1.1, 1.1, -1.1, '|\Psi^+\rangle');
    text(1.1, -1.1, 1.1, '|\Phi^+\rangle');
    text(-1.1, 1.1, 1.1, '|\Phi^-\rangle');
    text(-1.2, -1.2, -1.2, '|\Psi^-\rangle');
    ind = [6, 11, 16];
  case 'pos'
    title('pos correlations')
    xlabel('ZX');
    ylabel('XY');
    zlabel('YZ');
    text(1.1, -1.1, 1.1, '|y+,0\rangle +|y-,1\rangle');
    text(-1.1, 1.1, 1.1, '|y+,0\rangle -|y-,1\rangle');
    text(1.1, 1.1, -1.1, '|y-,0\rangle +|y+,1\rangle');
    text(-1.2, -1.2, -1.2, '|y-,0\rangle -|y+,1\rangle');
    ind = [8, 10, 15];
  case 'neg'
    title('neg correlations')
    xlabel('XZ');
    ylabel('YX');
    zlabel('ZY');
    text(1.1, 1.1, -1.1, '|0,y-\rangle +|1,y+\rangle');
    text(-1.1, 1.1, 1.1, '|0,y+\rangle -|1,y-\rangle');
    text(1.1, -1.1, 1.1, '|0,y+\rangle +|1,y-\rangle');
    text(-1.2, -1.2, -1.2, '|0,y-\rangle -|1,y+\rangle');
    ind = [14, 7, 12];
  case 'none'
    ind = [];
  otherwise
    error('Unknown set of correlations.')
end
