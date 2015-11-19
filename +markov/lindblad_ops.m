function [A, H_LS] = lindblad_ops(H, D, baths)
% LINDBLAD_OPS  Lindblad operators for a Born-Markov master equation.
%  [A, H_LS] = lindblad_ops(H, D, B)
%
%  Builds the Lindblad operators corresponding to a
%  base Hamiltonian H and a (Hermitian) interaction operator D
%  coupling the system to bath B.
%
%  Returns A == {A_i * sqrt(TU)}_i and H_LS * TU / \hbar,
%  where A_i are the Lindblad operators and H_LS is the Lamb shift.
%
%  B can also be a cell vector of baths, in which case D has to be
%  a cell vector of the corresponding interaction operators.

% Ville Bergholm 2009-2015


global qit

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
[dH, V] = markov.ops(H, D);

H_LS = 0;

for n=1:n_baths
  b = baths{n};
  ind = 1;
  for k=1:length(dH)
    % first the positive energy shift
    [g, s] = b.corr(dH(k));
    % is the dissipation significant?
    if abs(g) >= qit.tol
        A{n, ind} = sqrt(g) * V{n,k};
        %NA(n, ind) = norm(A{n, ind}, 'fro'); % how significant is this op?
        ind = ind+1;
    end
    % contribution to Lamb shift
    H_LS = H_LS +s * V{n,k}' * V{n,k};

    if dH(k) == 0
        % no negative shift
        continue
    end

    % now the corresponding negative energy shift
    [g, s] = b.corr(-dH(k));
    if abs(g) >= qit.tol
        A{n, ind} = sqrt(g) * V{n,k}';   % note the difference here, V(-omega) = V'(omega)
        %NA(n, ind) = norm(A{n, ind}, 'fro');
        ind = ind+1;
    end
    H_LS = H_LS +s * V{n,k} * V{n,k}'; % here too
  end
end
%NA


% TODO ops for different baths can be combined into a single basis,
% N^2-1 ops max in total
