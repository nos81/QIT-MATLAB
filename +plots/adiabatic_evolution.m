function adiabatic_evolution(t, st, H_func)
% ADIABATIC_EVOLUTION  Adiabatic evolution plot.
%  adiabatic_evolution(t, st, H_func)
%
%  Input: vector t of time instances, cell vector st of states corresponding
%  to the times and time-dependant Hamiltonian function handle H_func.
%
%  Plots the energies of the eigenstates of H_func(t(k)) as a function of t(k),
%  and the overlap of st{k} with the n lowest final Hamiltonian eigenstates. 
%  Useful for illustrating adiabatic evolution.

% Jacob D. Biamonte 2008
% Ville Bergholm 2009


n = 4;

T = t(end);
H = H_func(T);

% find the n lowest eigenstates of the final Hamiltonian
[v,d] = eig(H);
[S,I] = sort(diag(d), 'ascend');
for j=1:n
  lowest{j} = state(v(:, I(j)));
end

for k=1:length(t)
  tt = t(k);
  H = H_func(tt);
  energies(:,k) = sort(real(eig(H)), 'ascend');

  for j=1:n
    overlaps(j,k) = fidelity(lowest{j}, st{k})^2; % squared overlap with lowest final states
  end
end


subplot(2,1,1);
plot(t/T, energies);
grid on;
%set(gca, 'XTick', [0:0.25*t:t]);
%set(gca, 'XTickLabel', {'s = 0';'s = 0.25';'s = 0.5';'s = 0.75'; 's = 1'});
%set(gca, 'YTick', [min(min(energies)):0.25*m:max(max(energies))])
%set(gca,'YTickLabel',{'0';'0.5';'1';'1.5'; '2'})
title('Energy spectrum');
xlabel('Adiabatic time');
ylabel('Energy');
axis([0, 1, min(min(energies)), max(max(energies))]);


subplot(2,1,2);
plot(t/T, overlaps); %, 'LineWidth', 1.7);
grid on;
%set(gca, 'XTick', [0:0.25*t:t]);
%set(gca, 'XTickLabel', {'s = 0';'s = 0.25';'s = 0.5';'s = 0.75'; 's = 1'});
%set(gca, 'YTick', [0:0.5:1.1]);
%set(gca, 'YTickLabel', {'0';'0.5'; '1'});
title('Squared overlap of current state and final eigenstates');
xlabel('Adiabatic time');
ylabel('Probability');
%ylabel('|<l(s)|\psi(s)>|^2');
legend('|0\rangle','|1\rangle','|2\rangle','|3\rangle')
%axis([0, 1, 0, 1]);
axis([0, 1, 0, max(max(overlaps))]);
