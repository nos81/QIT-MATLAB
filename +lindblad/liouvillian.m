function L = liouvillian(H, D, baths)
% LINDBLAD/LIOUVILLIAN  Liouvillian superoperator for a Born-Markov system.
%
%  L = liouvillian(H, D, b)
%
%  Builds the Liouvillian superoperator L corresponding to a base Hamiltonian H
%  and a (hermitian) interaction operator D coupling the system to bath b.
%
%  Returns L/omega0, which includes the system Hamiltonian, the Lamb shift,
%  and the Lindblad dissipator.
%
%  b can also be a cell vector of baths, in which case D has to be
%  a cell vector of the corresponding interaction operators.

% Ville Bergholm 2009


% lindblad ops
[dH, X] = lindblad.ops(H, D);


if (~iscell(baths))
  baths = {baths}; % needs to be a cell array, even if it has just one element
end
n_baths = length(baths); % number of baths
% TODO the baths must share the same omega0!

iH_LS = 0;
acomm = 0;
diss = 0;

for n=1:n_baths
  A = X(n,:); % ops for the n'th bath
  b = baths{n};

  % we build the liouvillian in a funny order to be a bit more efficient

  % dH == 0 terms
  [g, s] = corr(b, 0);
  temp = A{1}' * A{1};

  iH_LS = iH_LS +i * s * temp; % Lamb shift
  acomm = acomm -0.5 * g * temp; % anticommutator
  diss = diss +lrmul(g * A{1}, A{1}'); % dissipator (part)

  for k=2:length(dH)
    % first the positive energy shift
    [g, s] = corr(b, dH(k));
    temp = A{k}' * A{k};

    iH_LS = iH_LS +i * s * temp;
    acomm = acomm -0.5 * g * temp;
    diss = diss +lrmul(g * A{k}, A{k}');

    % now the corresponding negative energy shift
    [g, s] = corr(b, -dH(k));
    temp = A{k} * A{k}'; % note the difference here, A(-omega) = A'(omega)

    iH_LS = iH_LS +i * s * temp;
    acomm = acomm -0.5 * g * temp;
    diss = diss +lrmul(g * A{k}', A{k}); % here too
  end
end

iH_LS = iH_LS +i*H; % include the system Hamiltonian

L = lmul(acomm -iH_LS) +rmul(acomm +iH_LS) +diss;
