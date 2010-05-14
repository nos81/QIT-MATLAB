function h = plot_makhlin_2q(sdiv, tdiv)
% PLOT_MAKHLIN_2Q  Plots the set of two-qubit gates in the space of Makhlin invariants.
%  h = plot_makhlin_2q(sdiv, tdiv)
%
%  Plots the set of two-qubit gates in the space of Makhlin
%  invariants, returns the surface handle h.
%
%  The input parameters are the s and t divisions of the mesh.

% Ville Bergholm 2006-2010


if (nargin < 2)
  tdiv = 41;
  if (nargin < 1)
    sdiv = 41;
  end
end

s = linspace(0, pi,   sdiv);
t = linspace(0, pi/2, tdiv);

% more efficient than meshgrid
%g1 = kron(cos(s).^2, cos(t).^4) - kron(sin(s).^2, sin(t).^4);
%g2 = 0.25*kron(sin(2*s), sin(2*t).^2);
%g3 = 4*g1 - kron(cos(2*s), cos(2*t).^2);
%S = kron(s, ones(size(t)));
%T = kron(ones(size(s)), t);

% canonical coordinate plane (s, t, t) gives the entire surface of the set of gate equivalence classes
[S, T] = meshgrid(s, t);
c = [S(:), T(:), T(:)];
G = invariant.makhlin(c);
G = reshape(G, sdiv, tdiv, 3);
C = invariant.max_concurrence(c);
C = reshape(C, sdiv, tdiv);

h = surf(G(:,:,1), G(:,:,2), G(:,:,3), C.^2);
%h = mesh(G(:,:,1), G(:,:,2), G(:,:,3), C.^2);
%h = waterfall(G(:,:,1), G(:,:,2), G(:,:,3), C.^2);
axis([-1 1 -0.5 0.5 -3 3]);
xlabel('g_1');
ylabel('g_2');
zlabel('g_3');
title('Makhlin stingray');

shading interp;
%hidden off;
set(gca, 'CLim', [0 1]); % color limits
colorbar

text(1.05, 0, 2.7, 'I');
text(-1.05, 0, -2.7, 'SWAP');
text(-0.1, 0, 1.2, 'CNOT');
text(0.1, 0, -1.2, 'DCNOT');
%text(0.1, 0.26, 0, 'SWAP^{1/2}');
%text(0, -0.26, 0, 'SWAP^{-1/2}');
