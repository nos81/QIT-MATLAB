function plot_weyl_2q()
% PLOT_WEYL_2Q  Plots the two-qubit Weyl chamber.
%  plot_weyl_2q()
%
%  Plots the Weyl chamber for the local invariants
%  of 2q gates.

% Ville Bergholm 2005-2010


hold off;
surf([0 0.5 1; 0 0.5 1], [0 0 0; 0 0.5 0], [0 0 0; 0 0.5 0]); alpha(0.2);
hold on;
surf([0 0.5; 0 0.5], [0 0.5; 0 0.5], [0 0; 0 0.5]); alpha(0.2);
surf([0.5 1; 0.5 1], [0.5 0; 0.5 0], [0 0; 0.5 0]); alpha(0.2);
%axis([0 1 0 0.5 0 0.5]);
axis equal;
xlabel('c_1/\pi');
ylabel('c_2/\pi');
zlabel('c_3/\pi');
title('Two-qubit Weyl chamber');

text(-0.05, -0.05, 0, 'I');
text(0.97, -0.05, 0, 'I');
text(0.45, 0.5, 0.55, 'SWAP');
text(0.45, -0.05, 0, 'CNOT');
text(0.45, 0.55, 0, 'DCNOT');
%text(0.20, 0.25, 0, 'SWAP^{1/2}');
%text(0.75, 0.25, 0, 'SWAP^{-1/2}');
