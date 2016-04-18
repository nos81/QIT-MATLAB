function nmr_sequences(seqs, titles)
% NMR_SEQUENCES  NMR robust control sequences demo.
%  nmr_sequences([seqs, titles])
%
%  Compares the performance of different single-qubit NMR control
%  sequences in the presence of systematic errors.
%  Plots the fidelity of each control sequence as a function
%  of both off-resonance error f and fractional pulse length error g.

% Ville Bergholm 2006-2016


fprintf('\n\n=== NMR control sequences for correcting systematic errors ===\n')

global qit;

if nargin < 1
    th = pi;
    seqs = {seq.nmr([th, 0]), seq.corpse(th), seq.scrofulous(th), seq.bb1(th), seq.knill(th)};
    titles = {'Plain \pi pulse', 'CORPSE', 'SCROFULOUS', 'BB1', 'Knill'};
else
    if nargin < 2
        titles = {'User-given seq'};
    end
end

% Pulse length/timing errors also affect the drift term, pulse strength errors don't.
strength_error = 1;
if strength_error
    error_type = 'pulse strength error';
else
    error_type = 'pulse length error';
end


% The two systematic error types we are interested here can be
% incorporated into the control sequence.

% off-resonance Hamiltonian (constant \sigma_z drift term) times -1i
offres_A = -1i * 0.5 * qit.sz;

psi = state('0'); % initial state

ns = length(seqs);
nf = 81;
ng = 71;
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
  if strength_error
      s_error.control = s.control * 1.1; % pulse strength error
  else
      s_error.tau = s.tau * 1.1; % pulse length error
  end

  % apply sequence on state psi, plot the evolution
  [out, t] = seq.propagate(psi, s_error, @bloch_vector);

  subplot(2,2,3);
  plot_state_trajectory(out);
  title([titles{q} ' evolution, ', error_type]);

  %==================================================
  s_error = s;

  for k=1:nf
      s_error.A = f(k) * offres_A; % off-resonance error
      for j=1:ng
          temp = 1 + g(j);
          if strength_error
              s_error.control = s.control * temp; % pulse strength error
          else
              s_error.tau = s.tau * temp; % pulse length error
          end
          fid(q, j, k) = u_fidelity(U, seq.seq2prop(s_error));
      end
  end

  subplot(2,2,[2 4]);
  [X,Y] = meshgrid(f,g);
  contour(X, Y, 1-squeeze(fid(q,:,:)));
  xlabel('Off-resonance error');
  ylabel(error_type);
  title([titles{q} ' fidelity']);
end

figure
plot(f, squeeze(fid(:, (ng+1)/2, :)));
xlabel('Off-resonance error');
ylabel('fidelity')
legend(titles)
grid on
axis([-1,1,0,1])

figure
plot(g, squeeze(fid(:, :, (nf+1)/2)));
xlabel(error_type);
ylabel('fidelity')
legend(titles)
grid on
axis([-1,1,0,1])
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
