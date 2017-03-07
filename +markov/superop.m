function L = superop(H, D, baths)
% SUPEROP  Liouvillian superoperator for a Born-Markov system.
%
%  L = superop(H, D, B)
%
%  Builds the Liouvillian superoperator L corresponding to a
%  system Hamiltonian H and a (Hermitian) interaction operator D
%  coupling the system to bath B in the weak coupling limit.
%
%  Returns L * TU, which includes the system Hamiltonian, the Lamb shift,
%  and the Lindblad dissipator.
%
%  B can also be a cell vector of baths, in which case D has to be
%  a cell vector of the corresponding interaction operators.

% Ville Bergholm 2009-2017


if ~iscell(baths)
  baths = {baths}; % needs to be a cell array, even if it has just one element
end
n_baths = length(baths); % number of baths

% make sure the baths have the same TU!
temp = baths{1}.TU;
for k=2:n_baths
  if baths{k}.TU ~= temp
    error('All the baths must have the same time unit!')
  end
end

% jump ops
[dH, X] = markov.ops(H, D);

iH_LS = 0;
acomm = 0;
diss = 0;

for n=1:n_baths
  A = X(n,:); % ops for the n'th bath
  b = baths{n};

  % we build the Liouvillian in a funny order to be a bit more efficient
  for k=1:length(dH)
    % first the positive energy shift
    [g, s] = b.corr(dH(k));
    temp = A{k}' * A{k};

    iH_LS = iH_LS +1i * s * temp;
    acomm = acomm -0.5 * g * temp;
    diss = diss +lrmul(g * A{k}, A{k}');

    if dH(k) == 0
        % no negative shift
        continue
    end

    % now the corresponding negative energy shift
    [g, s] = b.corr(-dH(k));
    temp = A{k} * A{k}'; % note the difference here, A(-omega) = A'(omega)

    iH_LS = iH_LS +1i * s * temp;
    acomm = acomm -0.5 * g * temp;
    diss = diss +lrmul(g * A{k}', A{k}); % here too
  end
end

iH_LS = iH_LS +1i*H; % include the system Hamiltonian

L = lmul(acomm -iH_LS) +rmul(acomm +iH_LS) +diss;
