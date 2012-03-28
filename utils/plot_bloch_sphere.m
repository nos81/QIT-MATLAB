function h = plot_bloch_sphere()
% PLOT_BLOCH_SPHERE  Bloch sphere plot.
%  h = plot_bloch_sphere()
%
%  Plots a Bloch sphere, a geometrical representation of the state space of a single qubit.
%  Pure states are on the surface of the sphere, nonpure states inside it.
%  The states |0> and |1> lie on the north and south poles of the sphere, respectively.

% Ville Bergholm  2005-2012
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
plot3(0,0,-1,'r.');
text(0, 0,  1.2, '|0\rangle');
text(0, 0, -1.2, '|1\rangle');

% equator
phi = linspace(0, 2*pi, 40);
plot3(cos(phi), sin(phi), zeros(size(phi)), 'k-');
end
