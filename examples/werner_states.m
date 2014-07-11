function res = werner_states(d)
% WERNER_STATES  Werner and isotropic states demo.
%
%  res = werner_states(d)
%

% Ville Bergholm 2014


fprintf('\n\n=== Werner and isotropic states ===\n\n')

if nargin < 1
    d = 2;
end

% cover both Werner ([0,1]) and the dual isotropic states
p = linspace(0, (d+1)/2, 200);
res = [];
for k=1:length(p)
    w = state(gate.werner(p(k), d));
    % corresponding isotropic state
    iso = w.ptranspose(1);
    res(k,:) = [w.purity(),   w.lognegativity(1),...
                iso.purity(), iso.lognegativity(1)];
end

figure();
plot(p, res);
hold on
% fully depolarized state
plot((d+1)/(2*d), 0, 'ko')
title(sprintf('Werner and isotropic states in d = %d', d))
xlabel('Werner state p');
legend('Werner purity', 'Werner lognegativity', 'Isotropic purity', 'Isotropic lognegativity');
grid on
