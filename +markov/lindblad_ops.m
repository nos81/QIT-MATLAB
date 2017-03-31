function [A, H_LS] = lindblad_ops(H, D, baths)
% LINDBLAD_OPS  Lindblad operators for a Born-Markov master equation.
%  [A, H_LS] = lindblad_ops(H, D, B)
%
%  Builds the Lindblad operators corresponding to a
%  system Hamiltonian H and a (Hermitian) interaction operator D
%  coupling the system to bath B in the weak coupling limit.
%
%  Returns A == {A_i * sqrt(TU)}_i and H_LS * TU / \hbar,
%  where A_i are the Lindblad operators and H_LS is the Lamb shift.
%
%  B can also be a cell vector of baths, in which case D has to be
%  a cell vector of the corresponding interaction operators.

% Ville Bergholm 2009-2017


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
[dH, V] = rotating_frame(H, D);

% The dH differences must not be too low for RWA to hold.
% The smallest dH difference determines \tau_S, the intrinsic time scale of
% the system, as defined in Breuer&Petruccione, chapter 3.3.1.
% For the RWA to work, \tau_S must be much shorter than \tau_R, the
% relaxation time of the system.

% tolerance for final transition frequency differences
tol_dH_warn = 1e-3;
for k=2:length(dH)
    if abs(dH(k)-dH(k-1)) < tol_dH_warn
        fprintf('Warning: The small difference between dH(%d) and dH(%d) may break the RWA.\n', k-1, k);
    end
end


H_LS = 0;

ind = 1;
for n=1:n_baths
  b = baths{n};
  for k=1:length(dH)
    % first the positive energy shift
    [g, s] = b.corr(dH(k));
    % is the dissipation significant?
    if abs(g) >= qit.tol
        A{ind} = sqrt(g) * V{k,n};
        %NA(ind) = norm(A{ind}, 'fro'); % how significant is this op?
        ind = ind+1;
    end
    % contribution to Lamb shift
    H_LS = H_LS +s * V{k,n}' * V{k,n};

    if dH(k) == 0
        % no negative shift
        continue
    end

    % now the corresponding negative energy shift
    [g, s] = b.corr(-dH(k));
    if abs(g) >= qit.tol
        A{ind} = sqrt(g) * V{k,n}';   % note the difference here, V(-omega) = V'(omega)
        %NA(ind) = norm(A{ind}, 'fro');
        ind = ind+1;
    end
    H_LS = H_LS +s * V{k,n} * V{k,n}'; % here too
  end
end
%NA


% TODO ops for different baths can be combined into a single basis,
% N^2-1 ops max in total
