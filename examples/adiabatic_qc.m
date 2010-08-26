function adiabatic_qc(H0, H1, s0, tmax)
% ADIABATIC_QC  Adiabatic quantum computing demo.
%  adiabatic_qc(H0, H1, s0, tmax)
%
%  This is a helper function for simulating the adiabatic quantum
%  algorithm of Farhi et al. and plotting the results.

%! E. Farhi et al., "Quantum Computation by Adiabatic Evolution", arXiv.org:quant-ph/0001106.
% Ville Bergholm 2009-2010


if (nargin < 4)
  tmax = 50; % how long the passage takes
end

H1_full = diag(H1); % into a full matrix

% adiabatic simulation
steps = tmax*10;
t = linspace(0, tmax, steps);

% linear path
H_func = @(t) (1-t/tmax)*H0 +(t/tmax)*H1_full;
res = propagate(s0, H_func, t);


% plots
% final state probabilities
figure;
plot_adiabatic_evolution(t, res, H_func);

figure;
plot(res{end});
title('Final state');

disp('Final Hamiltonian (diagonal):')
H1

disp('Measured result:')
[dummy, dummy, res] = measure(res{end});
display(res)
if (H1(find(res.data)) == 0)
  disp('Which is a valid solution!')
else
  disp('Which is not a solution!')
end
end
