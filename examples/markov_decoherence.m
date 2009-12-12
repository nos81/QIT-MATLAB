function [H, D, b] = markov_decoherence(T1, T2, b)
% MARKOV_DECOHERENCE  Decoherence demo.
%
%  [H, D, b] = markov_decoherence(T1, T2)
%
%  Given decoherence times T1 and T2, creates a markovian bath b
%  and a coupling operator D which reproduce them on a single-qubit system.

% Ville Bergholm 2009


global qit;

fprintf('\n\n=== Markovian decoherence ===\n')

if (nargin < 2)
  error('Both T1 and T2 must be given.')
end

omega0 = 1e9 % Hz
T = 3 % K

delta = 2.5 + rand; % qubit energy splitting

% setup the bath
if (nargin < 3)
  b = lindblad.bath('ohmic', omega0, T); % defaults
end

[b, H, D] = fit(b, delta, T1, T2);
L = lindblad.liouvillian(H, D, b);

t = 0:0.1:10;

% T1 demo
eq = 1/(1+exp(delta*b.scale)); % equilibrium rho_11

s = state('1', [2]); % qubit in the |1> state
out = propagate(s, L, t, @(x,h) ev(x, qit.p1));
figure
plot(t, cell2mat(out), 'r-', t, eq +(1-eq)*exp(-t/(T1*omega0)), 'b-.', 'LineWidth', 2);
xlabel('t \omega_0');
ylabel('probability');
axis([0 t(end) 0 1])
title('T_1: decay');
legend('P_1', 'exp(-t/T_1)')


% T2 demo
s = state('0', [2]);
s = u_propagate(s, R_y(pi/2)); % rotate to (|0>+|1>)/sqrt(2)
out = propagate(s, L, t, @(x,h) ev(u_propagate(x, R_y(-pi/2)), qit.p0));
figure
plot(t, cell2mat(out), 'r-', t, 0.5*(1+exp(-t/(T2*omega0))), 'b-', 'LineWidth', 2);
xlabel('t \omega_0');
ylabel('probability')
axis([0 t(end) 0 1])
title('T_2: dephasing');
legend('P_0', '0.5*(1+exp(-t/T_2))')
