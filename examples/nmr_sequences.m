function nmr_sequences
% NMR_SEQUENCES  NMR control sequences example.
%  nmr_sequences
%
%  Compares the performance of different single-qubit NMR control
%  sequences in the presence of systematic errors.
%  Plots the fidelity of each control sequence as a function
%  of both off-resonance error f and fractional pulse lenght error g.

%! Cummins et al., "Tackling systematic errors in quantum logic gates with composite rotations", PRA 67, 042308 (2003).
% Ville Bergholm 2006-2009 


fprintf('\n\n=== NMR control sequences for correcting systematic errors ===\n')

seqs = {[1 0 0 pi], seq.corpse(pi), seq.scrofulous(pi), seq.bb1(pi)};
titles = {'Plain \pi pulse', 'CORPSE', 'SCROFULOUS', 'BB1'};

psi = state('0'); % initial state

f = -1:0.05:1;
g = -1:0.08:1;
nf = length(f);
ng = length(g);

for q=1:length(seqs)
figure;

s = seqs{q};
U = seq.seq2prop(s); % target propagator

% in this simple example the errors can be fully included in the control sequence
%==================================================
s_error = s;
s_error(:,3) = s(:,3) +0.1; % off-resonance error

% apply sequence on state psi, plot the evolution
[out, t] = seq_propagate(psi, s_error, @bloch_vector);
a = cell2mat(out);
n = size(a, 2);

subplot(2,2,1);
plot_bloch_sphere();
plot3(a(1,:),a(2,:),a(3,:));
plot3(a(1,n),a(2,n),a(3,n), 'k.');
title([titles{q} ' evolution, off-resonance error']);

%==================================================
s_error = s;
s_error(:,end) = s(:,end)*1.1; % pulse lenght error

% apply sequence on state psi, plot the evolution
[out, t] = seq_propagate(psi, s_error, @bloch_vector);
a = cell2mat(out);
n = size(a, 2);

subplot(2,2,3);
plot_bloch_sphere();
plot3(a(1,:),a(2,:),a(3,:));
plot3(a(1,n),a(2,n),a(3,n), 'k.');
title([titles{q} ' evolution, pulse length error']);

%==================================================
s_error = s;
fid = [];

for k=1:nf
  s_error(:,3) = s(:,3) + f(k); % off-resonance error (constant \sigma_z interaction)
  for j=1:ng
    s_error(:,end) = s(:,end)*(1 + g(j)); % proportional pulse length error

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
