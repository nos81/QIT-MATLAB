function h = bloch_sphere(s)
% BLOCH_SPHERE  Bloch sphere plot.
%  h = bloch_sphere(s)
%
%  Plots a Bloch sphere, a geometrical representation of the state space of a single qubit.
%  Pure states are on the surface of the sphere, nonpure states inside it.
%  The states |0> and |1> lie on the north and south poles of the sphere, respectively.
%
%  s is a two dimensional state to be plotted.

% Ville Bergholm  2005-2010
% James Whitfield 2010


[X,Y,Z] = sphere(40);

hold off;
h = surf(X,Y,Z, 2*ones(41,41));
hold on;
shading flat
alpha(0.2)
axis square
xlabel('x');
ylabel('y');
zlabel('z');
plot3(0,0,1,'r.');
plot3(0,0,-1,'b.');

text(0, 0,  1.2, '|0\rangle');
text(0, 0, -1.2, '|1\rangle');
if nargin==1
    v = bloch_vector(s);
    quiver3(0, 0, 0, v(1), v(2), v(3), 0);
end