function h = wigner_function(W, a, b)
% WIGNER_FUNCTION  Wigner function plot.
%  h = wigner_function(W, x, y)
%
%  Plots the Wigner function given in the matrix W.
%  The vectors x and y define the coordinate grid.
%
% Ville Bergholm 2010


h = pcolor(a, b, W);
axis equal tight;
shading interp;
set(gca, 'CLim', [-1 1]);
colorbar;
colormap(asongoficeandfire(256));
xlabel('Re(\alpha)')
ylabel('Im(\alpha)')
title('Wigner function W(\alpha)')
