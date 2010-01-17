function [] = qubit_and_resonator()
% QUBIT_AND_RESONATOR  Simulates a qubit coupled to a resonator.
%
%  [] = qubit_and_resonator()
%
%  Simulates a qubit coupled to a resonator, replicating the experiment in the reference.

%! M. Hofheinz et al., "Synthesizing arbitrary quantum states in a superconducting resonator", Nature 459, 546-549 (2009), doi:10.1038/nature08005
% Ville Bergholm 2010


global qit


d_r = 5; % resonator dim

%b = lindblad.bath('ohmic', omega0, T);
omega0 = 1e9; % GHz

Omega     = 2*pi* 463e-3/25; %19e-3 % GHz
omega_r   = 2*pi* 6.570; % GHz
Delta_off = -25*Omega; % Hz

% Delta = omega_q -omega_r

%(omega_r -omega_off)^2/Omega^2

T1_q = 650e-9; % s
T2_q = 150e-9; % s

T1_r = 3.5e-6; % s
T2_r = 2*T1_r; % s

T = 0.025; % K

% qubit raising and lowering
sp = 0.5*(qit.sx -i*qit.sy);
sm = sp';
a = ho.ladder(d_r);


% rotating frame of the resonator
Hq = sp*sm;
Hint = Omega/2 * (kron(sp, a) +kron(sm, a'));
HOq = 0.5*(sp+sm);
HOr = 0.5*(a+a');
%Hr = omega_r * a'*a;

H = @(x) kron(x * Hq, eye(d_r)) +Hint;
% kron(Omega_q * HOq, eye(d_r))
% kron(qit.I, Omega_r * HOr)

s0 = state(0, [2, d_r]); % ground state
% s = tensor(state(0, [2]), ho.cstate(alpha, d_r)); % coherent state


% rabi test
out = [];
t = linspace(0, 300, 50);
ddd = linspace(0, 40, 50) * 2*pi*1e-3; % MHz
for k=1:length(ddd)
  %s = propagate(s, H0+microwave(omega_off), t);
  s = u_propagate(s0, kron(qit.sx, eye(d_r)));
  out(k,:) = cell2mat(propagate(s, H(ddd(k)), t, @(s,h) ev(s, kron(qit.p1, eye(d_r)))));
end
figure
pcolor(out);

figure
f = fft(out, [], 2);
pcolor(abs(fftshift(f, 2)));

return


% pumping
targ = zeros(d_r, 1);
targ(2) = 1;
%targ(3) = i;
targ = targ/norm(targ);
targ = tensor(state(0, 2), state(targ, d_r));

q = find(abs(targ.data), 1, 'last')-1; % highest excited level
prog = zeros(q, 3);

for k=q:-1:1
  prog(k, 1) = 1;
  prog(k, 2) = asin(abs(targ(k+1))/(sqrt(k)*Omega));
  s = propagate(s, H(0), prog(k, 2)); % S

  prog(k, 3) = 1;
end

s = s0;
for k=1:q
  % Q, S, Z
  s = propagate(s, H(Delta_off) +kron(Omega_q * HOq, eye(d_r)), t); % Q
  s = propagate(s, H(0), prog(k, 2)); % S
  s = propagate(s, H(Delta_off), t); Z
end

return




% readout
t = linspace(0, 10, 20);
out = propagate(s, H1, t, @(s,h) ev(s, kron(qit.p1, eye(d_r))));

