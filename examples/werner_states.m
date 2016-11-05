function werner_states(d)
% WERNER_STATES  Werner and isotropic states demo.
%  werner_states(d)
%
%  Plots some properties of d-dimensional family of Werner states and their
%  dual isotropic states as a function of the parameter p.

% Ville Bergholm 2014-2016


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
    res(k,:) = [w.purity(),...
                w.lognegativity(1),...
                iso.lognegativity(1)];
    %res2(k,1) = iso.purity();
end
%norm(res(:,1)-res2)  % approx. zero


figure();
leg = {'Werner/isotropic purity', 'Werner lognegativity', 'Isotropic lognegativity',...
       'maximally mixed state', 'maximally entangled generalized Bell state'};
plot(p, res);
hold on
% fully depolarized state
p_mix = (d+1)/(2*d);
plot(p_mix, 0, 'ko')
% generalized Bell state
p_bell = (d+1)/2;
plot(p_bell, 0, 'rs')
if d == 2
    % singlet state
    p_singlet = 0;
    plot(p_singlet, 0, 'bs')
    leg{end+1} = 'singlet state';
end
title(sprintf('Werner and isotropic states in d = %d', d))
xlabel('Werner state parameter p');
legend(leg);
grid on
