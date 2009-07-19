function l = makhlin_gates(sdiv, tdiv)
% MAKHLIN_GATES  Plots the surface of the set of 2q gates in the space of Makhlin invariants.
%  y = makhlin_gates(sdiv, tdiv)  % parameters are the s and t divisions of the mesh
%

% Ville Bergholm 2006-2008


if (nargin < 2)
  tdiv = 20;
  if (nargin < 1)
    sdiv = 20;
  end
end

ds = pi/sdiv;
s = 0:ds:pi;

dt = (pi/2)/tdiv;
t = (0:dt:pi/2)';

% more efficient than meshgrid
%g1 = kron(cos(s).^2, cos(t).^4) - kron(sin(s).^2, sin(t).^4);
%g2 = 0.25*kron(sin(2*s), sin(2*t).^2);
%g3 = 4*g1 - kron(cos(2*s), cos(2*t).^2);
%S = kron(s, ones(size(t)));
%T = kron(ones(size(s)), t);

% canonical coordinate plane (s, t, t) gives the entire surface of the set of gate equivalence classes
[S, T] = meshgrid(s, t);
[G1, G2, G3] = invariant.makhlin(S, T, T);
C = max_concurrence(S, T, T);

l = surf(G1, G2, G3, C.^2);
%l = mesh(G1, G2, G3, C);
%l = waterfall(G1, G2, G3, C);
axis([-1 1 -0.5 0.5 -3 3]);
xlabel('g1');
ylabel('g2');
zlabel('g3');
title('Makhlin stingray');
text(1.05, 0, 2.7, 'I');
text(-1.05, 0, -2.7, 'SWAP');
text(-0.1, 0, 1.2, 'CNOT');
text(0.1, 0, -1.2, 'DCNOT');
%text(0.1, 0.26, 0, 'SWAP^{1/2}');
%text(0, -0.26, 0, 'SWAP^{-1/2}');
shading interp;
%hidden off;
