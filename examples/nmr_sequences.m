function nmr_sequences(seqs, titles)
% NMR_SEQUENCES  NMR control sequences example.
%  nmr_sequences([seqs, titles])
%
%  Compares the performance of different single-qubit NMR control
%  sequences in the presence of systematic errors.
%  Plots the fidelity of each control sequence as a function
%  of both off-resonance error f and fractional pulse lenght error g.

%! Cummins et al., "Tackling systematic errors in quantum logic gates with composite rotations", PRA 67, 042308 (2003).
% Ville Bergholm 2006-2011


fprintf('\n\n=== NMR control sequences for correcting systematic errors ===\n')

global qit;

if nargin < 1
    seqs = {seq.nmr([pi, 0]), seq.corpse(pi), seq.scrofulous(pi), seq.bb1(pi)};
    titles = {'Plain \pi pulse', 'CORPSE', 'SCROFULOUS', 'BB1'};
else
    if nargin < 2
        titles = {'User-given seq'};
    end
end

psi = state('0'); % initial state

f = -1:0.05:1;
g = -1:0.08:1;
nf = length(f);
ng = length(g);

for q=1:length(seqs)
figure;

s = seqs{q};
U = seq.seq2prop(s); % target propagator

% The two systematic error types we are interested here can be
% incorporated into the control sequence.

%==================================================
s_error = s;
s_error.A = 0.1 * 1i * 0.5 * qit.sz; % add off-resonance error (constant \sigma_z drift term)

% apply sequence on state psi, plot the evolution
[out, t] = seq.propagate(psi, s_error, @bloch_vector);

subplot(2,2,1);
plot_state_trajectory(out);
title([titles{q} ' evolution, off-resonance error']);

%==================================================
s_error = s;
s_error.tau = s.tau * 1.1; % proportional pulse lenght error

% apply sequence on state psi, plot the evolution
[out, t] = seq.propagate(psi, s_error, @bloch_vector);

subplot(2,2,3);
plot_state_trajectory(out);
title([titles{q} ' evolution, pulse length error']);

%==================================================
s_error = s;
fid = [];

for k=1:nf
  s_error.A = f(k) * 1i * 0.5 * qit.sz; % off-resonance error
  for j=1:ng
    s_error.tau = s.tau * (1 + g(j)); % pulse length error
    fid(j, k) = u_fidelity(U, seq.seq2prop(s_error));
  end
end

subplot(2,2,[2 4]);
[X,Y] = meshgrid(f,g);
contour(X,Y,1-fid);
%surf(X,Y,1-fid);
xlabel('Off-resonance error');
ylabel('Pulse length error');
title([titles{q} ' fidelity']);
end
end


function F = u_fidelity(a,b)
% fidelity of two unitary rotations, [0,1]

  F = 0.5*abs(trace(a'*b));
end
