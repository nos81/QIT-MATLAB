function nmr_sequences(seqs, titles)
% NMR_SEQUENCES  NMR control sequences example.
%  nmr_sequences([seqs, titles])
%
%  Compares the performance of different single-qubit NMR control
%  sequences in the presence of systematic errors.
%  Plots the fidelity of each control sequence as a function
%  of both off-resonance error f and fractional pulse lenght error g.

%! Cummins et al., "Tackling systematic errors in quantum logic gates with composite rotations", PRA 67, 042308 (2003).
% Ville Bergholm 2006-2015


fprintf('\n\n=== NMR control sequences for correcting systematic errors ===\n')

global qit;

if nargin < 1
    th = pi;
    seqs = {seq.nmr([th, 0]), seq.corpse(th), seq.scrofulous(th), seq.bb1(th)};
    titles = {'Plain \pi pulse', 'CORPSE', 'SCROFULOUS', 'BB1'};
else
    if nargin < 2
        titles = {'User-given seq'};
    end
end

% The two systematic error types we are interested here can be
% incorporated into the control sequence.

% off-resonance Hamiltonian (constant \sigma_z drift term) times -1i
offres_A = -1i * 0.5 * qit.sz;

psi = state('0'); % initial state

ns = length(seqs);
nf = 41;
ng = 31;
f = linspace(-1, 1, nf);
g = linspace(-1, 1, ng);

fid = zeros(ns, ng, nf);
for q=1:ns
  figure;

  s = seqs{q};
  U = seq.seq2prop(s); % target propagator

  %==================================================
  s_error = s;
  s_error.A = 0.1 * offres_A; % add off-resonance error

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

  for k=1:nf
      s_error.A = f(k) * offres_A; % off-resonance error
      for j=1:ng
          s_error.tau = s.tau * (1 + g(j)); % pulse length error
          fid(q, j, k) = u_fidelity(U, seq.seq2prop(s_error));
      end
  end

  subplot(2,2,[2 4]);
  [X,Y] = meshgrid(f,g);
  contour(X, Y, 1-squeeze(fid(q,:,:)));
  xlabel('Off-resonance error');
  ylabel('Pulse length error');
  title([titles{q} ' fidelity']);
end

figure
plot(f, squeeze(fid(:, (ng+1)/2, :)));
xlabel('Off-resonance error');
ylabel('fidelity')
legend(titles)
grid on

figure
plot(g, squeeze(fid(:, :, (nf+1)/2)));
xlabel('Pulse length error');
ylabel('fidelity')
legend(titles)
grid on
end


function F = u_fidelity(a,b)
% fidelity of two unitary rotations, [0,1]

  F = 0.5*abs(trace(a'*b));
end


function F = s_fidelity(a,b)
% fidelity of two state transfers, [0,1]
% essentially just compare the first columns of the propagators

  temp = a'*b;
  F = abs(temp(1,1));
end
