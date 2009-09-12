function L = liouvillian(H, D, b)
% LINDBLAD/LIOUVILLIAN  Liouvillian superoperator for a Born-Markov system.
%
%  L = liouvillian(H, D, b)
%
%  Builds the Liouvillian superoperator L corresponding to a base Hamiltonian H
%  and a (hermitian) interaction operator D coupling the system to bath b.
%
%  Returns L/omega0, which includes the system Hamiltonian, the Lamb shift,
%  and the Lindblad dissipator.

% Ville Bergholm 2009


% lindblad ops
[dH, A] = lindblad.ops(H, D);


% we build the liouvillian in a funny order to be a bit more efficient

% dH == 0 terms
[g, s] = lindblad.bath_corr(b, 0);
temp = A{1}' * A{1};

iH_LS = i * s * temp; % Lamb shift
acomm = -0.5 * g * temp; % anticommutator
diss = lrmul(g * A{1}, A{1}'); % dissipator (part)

for k=2:length(dH)
  % first the positive energy shift
  [g, s] = lindblad.bath_corr(b, dH(k));
  temp = A{k}' * A{k};

  iH_LS = iH_LS +i * s * temp;
  acomm = acomm -0.5 * g * temp;
  diss = diss +lrmul(g * A{k}, A{k}');

  % now the corresponding negative energy shift
  [g, s] = lindblad.bath_corr(b, -dH(k));
  temp = A{k} * A{k}'; % note the difference here, A(-omega) = A'(omega)

  iH_LS = iH_LS +i * s * temp;
  acomm = acomm -0.5 * g * temp;
  diss = diss +lrmul(g * A{k}', A{k}); % here too
end

iH_LS = iH_LS +i*H; % include the system Hamiltonian

L = lmul(acomm -iH_LS) +rmul(acomm +iH_LS) +diss;
