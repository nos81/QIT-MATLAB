function h = bloch_sphere(s)
% BLOCH_SPHERE  Bloch sphere plot.
%  h = bloch_sphere(s)
%
%  Plots a Bloch sphere, a geometrical representation of the state space of a single qubit.
%  Pure states are on the surface of the sphere, nonpure states inside it.
%  The states |0> and |1> lie on the north and south poles of the sphere, respectively.
%
%  s is a two dimensional state to be plotted.

% James Whitfield 2010
% Ville Bergholm  2005-2009


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

text(0, 0,  1.1, '|0>');
text(0, 0, -1.1, '|1>');
if nargin==1
    vectarrow([0;0;0],bloch_vector(s));
end
