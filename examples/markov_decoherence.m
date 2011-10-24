function [H, D, B] = markov_decoherence(T1, T2, B)
% MARKOV_DECOHERENCE  Decoherence demo.
%
%  [H, D, B] = markov_decoherence(T1, T2)
%
%  Given decoherence times T1 and T2, creates a markovian bath B
%  and a coupling operator D which reproduce them on a single-qubit system.

% Ville Bergholm 2009-2010


global qit;

fprintf('\n\n=== Markovian decoherence in a qubit ===\n')

if (nargin < 2)
  error('Both T1 and T2 must be given.')
end

omega0 = 2*pi* 1e9; % Hz
T = 1 % K

delta = 3 + 3*rand; % qubit energy splitting (GHz)

% setup the bath
if (nargin < 3)
  B = markov.bath('ohmic', omega0, T); % defaults
end

% find the correct qubit-bath coupling
[H, D] = fit(B, delta, T1*omega0, T2*omega0)

L = markov.superop(H, D, B);
t = linspace(0, 10, 200);


% T1 demo
eq = 1/(1+exp(delta*B.scale)) % equilibrium rho_11

s = state('1', [2]); % qubit in the |1> state
out = propagate(s, L, t, @(x,h) ev(x, qit.p1));
figure
subplot(2, 1, 1)
plot(t, cell2mat(out), 'r-', t, eq +(1-eq)*exp(-t/(T1*omega0)), 'b-.',...
     [0, t(end)], [eq, eq], 'k:', 'LineWidth', 2);
xlabel('t \omega_0');
ylabel('probability');
axis([0 t(end) 0 1])
title('T_1: relaxation');
legend('P_1', 'P_{1}^{eq} + (1-P_{1}^{eq}) exp(-t/T_1)')


% T2 demo
s = state('0', [2]);
s = u_propagate(s, R_y(pi/2)); % rotate to (|0>+|1>)/sqrt(2)
out = propagate(s, L, t, @(x,h) ev(u_propagate(x, R_y(-pi/2)), qit.p0));
subplot(2, 1, 2);
plot(t, cell2mat(out), 'r-', t, 0.5*(1+exp(-t/(T2*omega0))), 'b-.', 'LineWidth', 2);
xlabel('t \omega_0');
ylabel('probability')
axis([0 t(end) 0 1])
title('T_2: dephasing');
legend('P_0', '(1 + exp(-t/T_2)) / 2')
