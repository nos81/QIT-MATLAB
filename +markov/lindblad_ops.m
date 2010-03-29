function [L, H_LS] = lindblad_ops(H, D, baths)
% MARKOV/LINDBLAD_OPS
%  [L, H_LS] = lindblad_ops(H, D, B)
%
%  Builds the Lindblad operators corresponding to a
%  base Hamiltonian H and a (Hermitian) interaction operator D
%  coupling the system to bath B.
%
%  Returns L == {A_i/omega0}_i and H_LS/(\hbar * omega0),
%  where A_i are the Lindblad operators and H_LS is the Lamb shift.
%
%  B can also be a cell vector of baths, in which case D has to be
%  a cell vector of the corresponding interaction operators.

% Ville Bergholm 2009-2010


if (~iscell(baths))
  baths = {baths}; % needs to be a cell array, even if it has just one element
end
n_baths = length(baths); % number of baths

% make sure the baths have the same omega0!
temp = baths{1}.omega0;
for k=2:n_baths
  if (baths{k}.omega0 ~= temp)
    error('All the baths must have the same energy scale omega0!')
  end
end

% jump ops
[dH, A] = markov.ops(H, D);

H_LS = 0;

for n=1:n_baths
  b = baths{n};

  % dH == 0 terms
  [g, s] = corr(b, 0);
  L{n,1} = sqrt(g) * A{n,1};
  H_LS = H_LS +s * A{n,1}' * A{n,1}; % Lamb shift

  for k=2:length(dH)
    % first the positive energy shift
    [g, s] = corr(b, dH(k));
    L{n, 2*k-2} = sqrt(g) * A{n,k};
    H_LS = H_LS +s * A{n,k}' * A{n,k};

    % now the corresponding negative energy shift
    [g, s] = corr(b, -dH(k));
    L{n, 2*k-1} = sqrt(g) * A{n,k}';   % note the difference here, A(-omega) = A'(omega)
    H_LS = H_LS +s * A{n,k} * A{n,k}'; % here too
  end
end

% TODO ops for different baths can be combined into a single basis,
% N^2-1 ops max in total
