function plot_adiabatic_evolution(t, st, H_func, n)
% PLOT_ADIABATIC_EVOLUTION  Adiabatic evolution plot.
%  plot_adiabatic_evolution(t, st, H_func)
%
%  Input: vector t of time instances, cell vector st of states corresponding
%  to the times and time-dependant Hamiltonian function handle H_func.
%
%  Plots the energies of the eigenstates of H_func(t(k)) as a function of t(k),
%  and the overlap of st{k} with the n lowest final Hamiltonian eigenstates. 
%  Useful for illustrating adiabatic evolution.

% Jacob D. Biamonte 2008
% Ville Bergholm 2009-2010


T = t(end);
H = H_func(T);

if (nargin < 4)
  n = 4;
end

n = min(n, length(H));

% find the n lowest eigenstates of the final Hamiltonian
[v,d] = eigs(H, n, 'SA');
[S,I] = sort(diag(d), 'ascend');
for j=1:n
  lowest{j} = state(v(:, I(j)));
end
% TODO with degenerate states these are more or less random linear combinations of the basis states... overlaps are not meaningful

for k=1:length(t)
  tt = t(k);
  H = H_func(tt);
  energies(:,k) = sort(real(eig(full(H))), 'ascend');

  for j=1:n
    overlaps(j,k) = fidelity(lowest{j}, st{k})^2; % squared overlap with lowest final states
  end
end


subplot(2,1,1);
plot(t/T, energies);
grid on;
title('Energy spectrum');
xlabel('Adiabatic time');
ylabel('Energy');
axis([0, 1, min(min(energies)), max(max(energies))]);


subplot(2,1,2);
plot(t/T, overlaps); %, 'LineWidth', 1.7);
grid on;
title('Squared overlap of current state and final eigenstates');
xlabel('Adiabatic time');
ylabel('Probability');
temp = char([]);
for k=1:n
  temp(k,:) = sprintf('|%d\\rangle', k-1);
end
legend(temp);
axis([0, 1, 0, 1]);
%axis([0, 1, 0, max(max(overlaps))]);
