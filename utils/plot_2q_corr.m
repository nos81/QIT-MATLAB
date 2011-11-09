function h = plot_2q_corr(s)
% PLOT_2Q_CORR  Plots the diagonal correlations simplex of a two-qubit state.
%  h = plot_2q_corr([s])
%
%  Plots the geometrical representation of the set of allowed diagonal
%  correlations in a two-qubit state. The set is a tetrahedron,
%  the vertices of which correspond to the four Bell states.
%
%  s is an optional two-qubit state to be plotted in the simplex.

% Ville Bergholm  2011


% vertices and faces
v = [1, 1, -1;  1, -1, 1;  -1, 1, 1; -1, -1, -1];
f = [1,2,4; 2,3,4; 3,1,4; 3,2,1];

hold off;
% tetrahedron
h = patch('Faces',f, 'Vertices',v, 'FaceVertexCData', 0, 'FaceColor', 'flat');
hold on;
alpha(0.2)
axis equal
grid on
xlabel('xx');
ylabel('yy');
zlabel('zz');

% mark vertices
plot3(0, 0, 0, 'r.'); % center
plot3(v(:,1), v(:,2), v(:,3), 'r.');
text(1.1, 1.1, -1.1, '|\Psi^+\rangle');
text(1.1, -1.1, 1.1, '|\Phi^+\rangle');
text(-1.1, 1.1, 1.1, '|\Phi^-\rangle');
text(-1.1, -1.1, -1.1, '|\Psi^-\rangle');

if nargin == 1
    b = diag(bloch_vector(s, true));
    plot3(b(2), b(3), b(4), 'k.');
end
