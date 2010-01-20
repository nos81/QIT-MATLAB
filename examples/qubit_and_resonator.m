function qubit_and_resonator(d_r)
% QUBIT_AND_RESONATOR  Demo: Qubit coupled to a microwave resonator.
%
%  qubit_and_resonator()
%
%  Simulates a qubit coupled to a microwave resonator.
%  Reproduces plots from the experiment in the reference.

%! M. Hofheinz et al., "Synthesizing arbitrary quantum states in a superconducting resonator", Nature 459, 546-549 (2009), doi:10.1038/nature08005
% Ville Bergholm 2010


global qit

if (nargin < 1 || d_r < 10)
  d_r = 10; % resonator dim
end

sx = qit.sx;
sy = qit.sy;
p1 = qit.p1;

%b = lindblad.bath('ohmic', omega0, T); % TODO noise
omega0 = 1e9; % Hz

Omega     = 2*pi* 19e-3; % GHz
omega_r   = 2*pi* 6.570; % GHz
Delta_off = -2*pi* 463e-3; % GHz
Omega_q   = 2*pi*0.5; % GHz, stronger than anything else


T1_q = 650e-9; % s
T2_q = 150e-9; % s

T1_r = 3.5e-6; % s
T2_r = 2*T1_r; % s

T = 0.025; % K

% qubit raising and lowering ops
sp = 0.5*(sx -i*sy);
sm = sp';

% resonator annihilation op
a = ho.ladder(d_r);


Hq = kron(eye(d_r), sp*sm);
%Hr = kron(omega_r * a'*a, qit.I);

Hint = Omega/2 * (kron(a, sp) +kron(a', sm));
HOqr = kron(eye(d_r), Omega_q*0.5*sx);
HOqi = kron(eye(d_r), Omega_q*0.5*sy);
%HOr = kron(Omega_r*0.5*(a+a'), qit.I);

% Hamiltonian, rotating frame defined by H0 = omega_r*Hq +Hr
H = @(D,Oqr,Oqi) D*Hq +Hint +Oqr*HOqr +Oqi*HOqi;

% readout: qubit in excited state?
readout = @(s,h) ev(s, kron(eye(d_r), p1));

s0 = state(0, [d_r, 2]); % ground state


%=================================
% rabi test

out = [];
t = linspace(0, 500, 100);
ddd = linspace(0, 2*pi* 40e-3, 100); % MHz
for k=1:length(ddd)
  s = propagate(s0, H(Delta_off, 1, 0), (2/Omega_q)*pi/2); % Q
  %s = u_propagate(s0, kron(eye(d_r), sx));
  out(k,:) = cell2mat(propagate(s, H(ddd(k), 0, 0), t, readout));
end
figure
pcolor(t, ddd/(2*pi*1e-3), out);
shading interp;
%shading flat;
xlabel('Interaction time \tau (ns)')
ylabel('Detuning, \Delta/(2\pi) (MHz)')
title('One photon Rabi-swap oscillations between qubit and resonator.')

%figure
%f = fft(out, [], 2);
%pcolor(abs(fftshift(f, 2)));


%=================================
% readout plot (not corrected for limited visibility)

sq = state(0, 2);

s1 = tensor(normalize(state(1, d_r) +state(3, d_r)), sq);
s2 = tensor(normalize(state(1, d_r) +i*state(3, d_r)), sq);

t = linspace(0, 350, 200);
out = [];
out(1,:) = cell2mat(propagate(s1, H(0, 0, 0), t, readout));
out(2,:) = cell2mat(propagate(s2, H(0, 0, 0), t, readout));


figure
plot(t, out(1,:), 'r-')
hold on
plot(t, out(2,:), 'b:')
xlabel('Interaction time \tau (ns)')
ylabel('P_e')
title('Resonator readout through qubit.')


%=================================
% state preparation and Wigner spectroscopy

function prog = demolish_state(targ)
% convert a desired resonator state ket into a program

% Ideal H without interaction
A = @(D,Oqr,Oqi) D*Hq +Oqr*HOqr +Oqi*HOqi;

q = n-1; % highest excited level
prog = zeros(q, 4);

for k=q:-1:1
  % |g,k> to |e,k-1> 
  dd = targ.data;
  prog(k, 4) = (angle(dd(2*k)) -angle(dd(2*k+1)) -pi/2 +2*pi) / -Delta_off;
  targ = propagate(targ, A(-Delta_off, 0, 0), prog(k, 4)); % Z

  dd = targ.data;
  prog(k, 3) = (2/(sqrt(k)*Omega))*atan2(abs(dd(2*k+1)), abs(dd(2*k)));
  targ = propagate(targ, -Hint, prog(k, 3)); % S

  % |e,k-1> to |g,k-1>
  dd = targ.data;
  phi = angle(dd(2*k)) -angle(dd(2*k-1)) +pi/2;
  prog(k, 2) = phi;
  prog(k, 1) = (2/Omega_q)*atan2(abs(dd(2*k)), abs(dd(2*k-1)));
  targ = propagate(targ, A(0, -cos(phi), -sin(phi)), prog(k, 1)); % Q
end
end


% state prep in reverse
if (false)
  % "voodoo cat"
  targ = zeros(d_r, 1);
  for k = 1:3:d_r
    targ(k) = 2^(k-1)/sqrt(factorial(k-1));
  end
else
  targ = [1 0 0 exp(i*pi/4) 0 0 1].';
end

n = length(targ);
targ = tensor(state([targ; zeros(d_r-n, 1)], d_r), state(0, 2));
targ = normalize(targ);

% calculate the pulse sequence for constructing targ
prog = demolish_state(targ);

s = s0; % ground state
%s = tensor(ho.coherent_state(0.5, d_r), state(0, 2)); % coherent state

for k=1:q
  % Q, S, Z
  s = propagate(s, H(0, cos(prog(k, 2)), sin(prog(k, 2))), prog(k, 1)); % Q
  s = propagate(s, H(0, 0, 0), prog(k, 3)); % S
  s = propagate(s, H(Delta_off, 0, 0), prog(k, 4)); % Z
end
fprintf('Overlap of prepared state with target state: %g\n', overlap(s, targ));

s = ptrace(s, 2);

figure
[W, a, b] = ho.wigner(s, [80 80]);
pcolor(a, b, W);
shading interp;
%shading flat;
xlabel('Re(\alpha)')
ylabel('Im(\alpha)')
title('Wigner function W(\alpha)')
end
