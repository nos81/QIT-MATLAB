function qubit_and_resonator(d_r)
% QUBIT_AND_RESONATOR  Demo: Qubit coupled to a microwave resonator.
%
%  qubit_and_resonator()
%
%  Simulates a qubit coupled to a microwave resonator.
%  Reproduces plots from the experiment in the reference.

%! M. Hofheinz et al., "Synthesizing arbitrary quantum states in a superconducting resonator", Nature 459, 546-549 (2009), doi:10.1038/nature08005
% Ville Bergholm 2010


fprintf('\n\n=== Qubit coupled to a single-mode microwave resonator. ===\n\n')

global qit

if (nargin < 1 || d_r < 10)
  d_r = 30; % resonator dim
end

sx = qit.sx;
sy = qit.sy;
p1 = qit.p1;

omega0 = 1e9; % Energy scale/\hbar, Hz
T = 0.025; % K

% omega_r = 2*pi* 6.570; % GHz, resonator angular frequency
% qubit-resonator detuning Delta(t) = omega_q(t) -omega_r

Omega     =  2*pi* 19e-3;  % GHz, qubit-resonator coupling
Delta_off = -2*pi* 463e-3; % GHz, detuning at off-resonance point
Omega_q   =  2*pi* 0.5;    % GHz, qubit microwave drive amplitude, stronger than anything else

% decoherence times
T1_q = 650e-9; % s
T2_q = 150e-9; % s

T1_r = 3.5e-6; % s
T2_r = 2*T1_r; % s

% heat bath and couplings
%bq = markov.bath('ohmic', omega0, T);
%[H, Dq] = fit(bq, Delta_off, T1_q, T2_q);
%[H, Dr] = fit_ho(bq, ???, T1_r, T2_r???);
%D = kron(eye(d_r), D); % +kron(Dr, qit.I);


%=================================
% Hamiltonians etc.

% qubit raising and lowering ops
sp = 0.5*(sx -i*sy);
sm = sp';

% resonator annihilation op
a = ho.ladder(d_r);
%Q = ho.position(d_r);
%P = ho.momentum(d_r);
% resonator identity op
I_r = eye(d_r);

Hq = kron(I_r, sp*sm);
%Hr = kron(omega_r * a'*a, qit.I);

Hint = Omega/2 * (kron(a, sp) +kron(a', sm));
HOq = @(Oq_a, Oq_p) kron(I_r, Oq_a*0.5*Omega_q*(cos(Oq_p)*sx +sin(Oq_p)*sy));
%HOr = @(Or_a, Or_p) kron(Or_a*Omega_r/sqrt(2)*(cos(Or_p)*Q +sin(Or_p)*P), qit.I);

% Hamiltonian, rotating frame defined by H0 = omega_r*Hq +Hr
H = @(D, Oq_a, Oq_p) D*Hq +Hint +HOq(Oq_a, Oq_p);

% readout: qubit in excited state?
readout = @(s,h) ev(s, kron(I_r, p1));

s0 = state(0, [d_r, 2]); % ground state


%=================================
% Rabi test

out = [];
t = linspace(0, 500, 100);
ddd = linspace(0, 2*pi* 40e-3, 100); % MHz
%L = markov.superop(H(0, 1, 0), D, bq);
L = H(0, 1, 0);
for k=1:length(ddd)
  s = propagate(s0, L, (2/Omega_q)*pi/2); % Q
  %s = u_propagate(s0, kron(I_r, sx));

  %LL = markov.superop(H(ddd(k), 0, 0), D, bq);
  LL = H(ddd(k), 0, 0);
  out(k,:) = cell2mat(propagate(s, LL, t, readout));
end

figure
pcolor(t, ddd/(2*pi*1e-3), out);
colormap(asongoficeandfire(256));
colorbar('northoutside')
set(gca, 'CLim', [0 1]);
shading interp;
xlabel('Interaction time \tau (ns)')
ylabel('Detuning, \Delta/(2\pi) (MHz)')
title('One photon Rabi-swap oscillations between qubit and resonator, P_e')

%figure
%f = fft(out, [], 2);
%pcolor(abs(fftshift(f, 2)));


%=================================
% state preparation

function [prog, t] = demolish_state(targ)
% convert a desired (possibly truncated) resonator state ket into a program
% state preparation in reverse

% Ideal H without interaction
A = @(D, Oq_a, Oq_p) D*Hq +HOq(Oq_a, Oq_p);

% resonator ket into a full normalized qubit+resonator state
n = length(targ);
targ = normalize(state([targ; zeros(d_r-n, 1)], d_r));
targ = tensor(targ, state(0, 2));
t = targ;

n = n-1; % highest excited level in resonator
prog = zeros(n, 4);

for k=n:-1:1
  % |k,g> to |k-1,e> 
  dd = targ.data;
  prog(k, 4) = (angle(dd(2*k)) -angle(dd(2*k+1)) -pi/2 +2*pi) / -Delta_off;
  targ = propagate(targ, A(-Delta_off, 0, 0), prog(k, 4)); % Z

  dd = targ.data;
  prog(k, 3) = (2/(sqrt(k)*Omega))*atan2(abs(dd(2*k+1)), abs(dd(2*k)));
  targ = propagate(targ, -Hint, prog(k, 3)); % S

  % |k-1,e> to |k-1,g>
  dd = targ.data;
  phi = angle(dd(2*k)) -angle(dd(2*k-1)) +pi/2;
  prog(k, 2) = phi;
  prog(k, 1) = (2/Omega_q)*atan2(abs(dd(2*k)), abs(dd(2*k-1)));
  targ = propagate(targ, A(0, -1, phi), prog(k, 1)); % Q
end
end


function s = prepare_state(prog)
% prepare a state according to the program

s = s0; % start with ground state
%s = tensor(ho.coherent_state(0.5, d_r), state(0, 2)); % coherent state

for k=1:size(prog, 1)
  % Q, S, Z
  s = propagate(s, H(0, 1, prog(k, 2)), prog(k, 1)); % Q
  s = propagate(s, H(0, 0, 0), prog(k, 3)); % S
  s = propagate(s, H(Delta_off, 0, 0), prog(k, 4)); % Z
end
end


%=================================
% readout plot (not corrected for limited visibility)

[prog] = demolish_state([0 1 0 1].'); % |1> + |3>
s1 = prepare_state(prog);

[prog] = demolish_state([0 1 0 i].'); % |1> + i|3>
s2 = prepare_state(prog);

t = linspace(0, 350, 200);
out = [];
out(1,:) = cell2mat(propagate(s1, H(0, 0, 0), t, readout));
out(2,:) = cell2mat(propagate(s2, H(0, 0, 0), t, readout));

figure
plot(t, out(1,:), 'b-', t, out(2,:), 'r-')
xlabel('Interaction time \tau (ns)')
ylabel('P_e')
title('Resonator readout through qubit.')
legend('|1\rangle + |3\rangle', '|1\rangle + i|3\rangle')


%=================================
% Wigner spectroscopy

if (false)
  % "voodoo cat"
  targ = zeros(d_r, 1);
  for k = 0:3:(d_r - 1)
    targ(k+1) = 2^k/sqrt(factorial(k));
  end
else
  targ = [1 0 0 exp(i*pi*3/8) 0 0 1].';
end

% calculate the pulse sequence for constructing targ
[prog, t] = demolish_state(targ);
s = prepare_state(prog);

disp('Trying to prepare the state')
display(t, 'short')
fprintf('Fidelity of prepared state with target state: %g\n', fidelity(s, t));
fprintf('Time required for state preparation: %g ns\n', sum(sum(prog(:, [1 3 4]))));

s = ptrace(s, 2);

figure
[W, a, b] = ho.wigner(s, [80 80], [-2.5 2.5 -2.5 2.5]);
pcolor(a, b, W);
axis equal tight;
shading interp;
set(gca, 'CLim', [-1 1]);
colorbar;
colormap(asongoficeandfire(256));
xlabel('Re(\alpha)')
ylabel('Im(\alpha)')
title('Wigner function W(\alpha)')
end
