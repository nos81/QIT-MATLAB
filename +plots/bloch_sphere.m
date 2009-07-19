function bloch_sphere()
% BLOCH_SPHERE  Bloch sphere plot.
%  bloch_sphere()
%
%  Plots a Bloch sphere, a geometrical representation of the state space of a single qubit.
%  Pure states are on the surface of the sphere, nonpure states inside it.
%  The states |0> and |1> are marked with red and blue dots, respectively.

% Ville Bergholm 2005-2009


[X,Y,Z] = sphere(40);

hold off;
surf(X,Y,Z, 2*ones(41,41))
hold on;
shading flat
alpha(0.2)
axis square
xlabel('x');
ylabel('y');
zlabel('z');
%title('Alice')
plot3(0,0,1,'r.');
plot3(0,0,-1,'b.');

