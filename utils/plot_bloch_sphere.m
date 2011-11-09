function h = plot_bloch_sphere(s)
% PLOT_BLOCH_SPHERE  Bloch sphere plot.
%  h = plot_bloch_sphere(s)
%
%  Plots a Bloch sphere, a geometrical representation of the state space of a single qubit.
%  Pure states are on the surface of the sphere, nonpure states inside it.
%  The states |0> and |1> lie on the north and south poles of the sphere, respectively.
%
%  s is an optional two-dimensional state to be plotted in the sphere.

% Ville Bergholm  2005-2011
% James Whitfield 2010


[X,Y,Z] = sphere(40);

hold off;
h = surf(X,Y,Z, 2*ones(41,41));
hold on;
shading flat
alpha(0.2)
axis equal
grid on
xlabel('x');
ylabel('y');
zlabel('z');

% poles
plot3(0,0,1,'r.');
plot3(0,0,-1,'b.');
text(0, 0,  1.2, '|0\rangle');
text(0, 0, -1.2, '|1\rangle');

% equator
phi = linspace(0, 2*pi, 40);
plot3(cos(phi), sin(phi), zeros(size(phi)), 'k-');

if nargin==1
    v = bloch_vector(s);
    quiver3(0, 0, 0, v(2), v(3), v(4), 0);
end
