function h = plot_pcolor(W, a, b, clim)
% PLOT_PCOLOR  Easy pseudocolor plot.
%  h = plot_pcolor(W, x, y [, clim])
%
%  Plots the 2D function given in the matrix W.
%  The vectors x and y define the coordinate grid.
%  clim is an optional parameter for color limits (default is [0 1]).
%
%  Returns the handle to the surface object.
%
% Ville Bergholm 2010


if (nargin < 4)
  clim = [0 1];
end

h = pcolor(a, b, W);
axis equal tight;
shading interp;
set(gca, 'CLim', clim);
colorbar;
colormap(asongoficeandfire(256));
