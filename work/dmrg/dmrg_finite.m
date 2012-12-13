function [E, block] = dmrg_finite(H_func, n, m, sweeps)
% DMRG_FINITE  DMRG method for solving finite 1d-nn quantum systems.
%  [E0, B] = dmrg_finite(H_func, n, m, sweeps)
%
%  Returns the approximate ground state energy for a
%  one-dimensional chain of n sites with nearest-neighbor
%  couplings, with a Hamiltonian defined by H_func, using a
%  finite-system density matrix renormalization group method.
%
%  Also returns the DMRG block structure B which contains a
%  description of the state.
%
%  [H, dim] = H_func([first, last]) should return the partial
%  Hamiltonian and dimension vector for the sites first:last.
%  [S1, S2, h1, h2] = H_func(s) should return the coupling
%  and local operators for sites s and s+1. 
%
%  At each truncation step, at most m states are kept (possibly
%  fewer near the ends of the chain).
%  sweeps is the number of back-and-forth sweeps along the chain.

%! White, Steven R., "Density matrix formulation for quantum renormalization groups", PRL 69, 2863--2866 (1992), doi:10.1103/PhysRevLett.69.2863
%! G. De Chiara et al., "Density Matrix Renormalization Group for Dummies", JCTN 5, 1277 (2008), doi:10.1166/jctn.2008.011

% Ville Bergholm 2010


% TODO we treat block almost as a global variable here (bad form!)
% try to eliminate embedded functions.

if (nargin < 4)
  sweeps = 3;
  if (nargin < 3)
    m = 5;
    if (nargin < 2)
      error('H_func and n are required.');
    end
  end
end

block = cell(1, n-1);

% block k:   L(1:k)|R(k+1:n), S1|S2
% contents: HL, HR
% vector S1: acting on the rightmost site represented by the L block
% vector S2: acting on the leftmost site represented by the R block
% h1, h2 likewise

% precompute ops (if constant, this is a waste)
for k=1:n-1
  [block{k}.S1, block{k}.S2, block{k}.h1, block{k}.h2] = H_func(k);    
end

b = 1 % initial block size

% initial L block
block{b}.P = 0;
H = H_func([1 b]);
block{b}.OL = speye(length(H));
block{b} = projectL(block{b}, H);

% initial R block
block{n-b}.P = 0;
H = H_func([n-b+1 n]);
block{n-b}.OR = speye(length(H));
block{n-b} = projectR(block{n-b}, H);


E = NaN(1,n-1);

skip = 0 % skip this many end sites in the sweeps TODO automatic

% initialization
if (1)
  % L*|*R scheme (heuristic)
  % grow the chain (only needed if n > 4)
  temp = ceil(n/2)-1;
  for k=b+1:temp
    % build block Hamiltonians
    prev = k-1;
    [HL, Ldim] = Lstar(block{prev});
    [HR, Rdim] = starR(block{n-prev});

    % NOTE the ** coupling is undefined until the blocks touch so we just make it up!
    % The superblock Hamiltonian is formed using the left coupling only.
    block{k}.state = []; % no guess
    [block{k}, E(k)] = combine(block{k}, HL, HR, [Ldim Rdim]);
    block{k} = truncate(block{k}, m);
    
    % copy OR transform so we can projectR
    block{n-k}.OR = block{k}.OR;
    
    block{k}   = projectL(block{k},   HL);
    block{n-k} = projectR(block{n-k}, HR);
  end

  % rest of sweep to right
  range = temp+1:n-b-skip;
  block{range(1)}.state = []; % no guess
  E(end+1,:) = sweep(range);
  
else
  % L*|**** scheme
  ns = 4; % number of sites in the right block

  for k=b+1:n-b-1  %n-1
    % build block Hamiltonians
    % the L* block has k sites in it
    [HL, Ldim] = Lstar(block{k-1});
  
    sites = min(ns, n-k);
    [HR, temp] = H_func([k+1, k+sites]);
    Rdim = [temp(1), length(HR)/temp(1)];
  
    block{k}.state = []; % no guess
    [block{k}, E(k)] = combine(block{k}, HL, HR, [Ldim Rdim]);
    block{k} = truncate(block{k}, m);

    block{k} = projectL(block{k}, HL);
    % projectR is unnecessary here
  end
end


% sweeps
rangeL = fliplr(b+skip:n-b-skip-1);
rangeR = b+skip+1:n-b-skip;
if (length(rangeL) == 0)
  return
end

fprintf('Sweeping')
for s=1:sweeps
  % so much trouble just to avoid double diagonalization at ends!
  q = rangeR(end);
  block{q-1}.state = guess_stateL(block{q}, block{q-1});
  E(end+1,:) = sweep(rangeL); % left
  
  q = rangeL(end);
  block{q+1}.state = guess_stateR(block{q}, block{q+1});
  E(end+1,:) = sweep(rangeR); % right
  
  fprintf('.')
end
fprintf(' done.\n')


function res = sweep(sites)
% L*|*R sweeping (also can handle L*|*, *|*R and *|*)
% TODO how about L*|R ?

  res = NaN(1,n-1);

  % this is OK with only one site too, because then direction does not matter.
  sweep_right = (sites(1) <= sites(end));

  for q=sites
    % after this iteration, block q will be replaced

    % build block Hamiltonians for L(q-1)* and *R(q+1)
    if (q==1)
        % just *|, no L
        Ldim = [1, length(block{q}.h1)];
        HL = block{q}.h1;
    else
        [HL, Ldim] = Lstar(block{q-1});
    end

    if (q == n-1)
        % just |*, no R
        Rdim = [length(block{q}.h2), 1];
        HR = block{q}.h2;
    else
        [HR, Rdim] = starR(block{q+1});
    end

    [block{q}, res(q)] = combine(block{q}, HL, HR, [Ldim Rdim]);
    block{q} = truncate(block{q}, m);
    
    % are we at the endpoint?
    if (q == sites(end))
      % endpoint, both projections, no guessing
      block{q} = projectR(block{q}, HR);
      block{q} = projectL(block{q}, HL);
      return
    end

    % trailing projection, forward guess
    if (sweep_right)
      block{q} = projectL(block{q}, HL);
      block{q+1}.state = guess_stateR(block{q}, block{q+1});
    else
      block{q} = projectR(block{q}, HR);
      block{q-1}.state = guess_stateL(block{q}, block{q-1});
    end
  end
end

end


% White's state prediction
function s = guess_stateR(B, C)
% make a guess towards the right
  dim = B.state.dim;
  O1 = B.OL';
  O2 = C.OR;
  c = length(C.h2);
  dd = [size(O1,1), length(C.h1), c, size(O2,1)/c];
  s = state(kron(speye(size(O1, 1)*dim(3)), O2)*(kron(O1, speye(prod(dim(3:4))))*B.state.data), dd);
end
      
function s = guess_stateL(B, A)
  dim = B.state.dim;
  O1 = B.OR';
  O2 = A.OL;
  a = length(A.h1);
  dd = [size(O2,1)/a, a, length(A.h2), size(O1,1)];
  s = state(kron(O2, speye(dim(2)*size(O1, 1)))*(kron(speye(prod(dim(1:2))), O1)*B.state.data), dd);
end


function C = coupling(S1, S2)
% C = \sum_k S1_k * S2_k
  C = sparse(0);
  for k=1:length(S1)
    C = C +kron(S1{k}, S2{k});
  end
end

function [H, dim] = Lstar(B)
% L* block
  dim = [length(B.HL), length(B.h2)];
  H = kron(B.HL, speye(dim(2))) +coupling(B.S1p, B.S2) +kron(speye(dim(1)), B.h2);
end

function [H, dim] = starR(B)
% *R block
  dim = [length(B.h1), length(B.HR)];
  H = kron(speye(dim(1)), B.HR) +coupling(B.S1, B.S2p) +kron(B.h1, speye(dim(2)));
end

function [B, E0] = combine(B, H_L, H_R, dim)
% NOTE: H_L and H_R are untruncated ops (unlike B.HL and B.HR)

  % superblock Hamiltonian
  H = kron(H_L, speye(length(H_R))) +kron(speye(length(H_L)), H_R)...
      +mkron(speye(dim(1)), coupling(B.S1, B.S2), speye(dim(4)));

  H = 0.5*(H+H'); % eliminate rounding errors

  % find target state TODO other targets?
  opts = struct();
  if (isa(B.state, 'state'))
    opts.v0 = B.state.data;
  end
  [v,E0] = eigs(H, 1, 'SA', opts); % ground state
  B.state = state(v, dim); % store it

  %if (~isscalar(guess))
  %  guess_fidelity = fidelity(guess, B.state) % test
  %end
end

function B = truncate(B, m)
  % Schmidt decompose block state (4 components) into two halves
  [d, u, v] = schmidt(B.state, [1 2]);

  keep = min(m, length(d)); % how many states to keep?
    
  B.P = 1 - sum(d(1:keep).^2); % truncation error
  % keep most influential states
  B.OL = u(:, 1:keep);
  B.OR = v(:, 1:keep);    
end

function B = projectL(B, H)
  B.HL = B.OL'*H*B.OL; % project the block left Hamiltonian
  
  % project S1 into the new block
  % we don't need h1 because it's already included in H
  temp = speye(length(B.OL)/length(B.S1{1}));
  for k=1:length(B.S1)
    B.S1p{k} = B.OL'*kron(temp, B.S1{k})*B.OL;
  end
end

function B = projectR(B, H)
  B.HR = B.OR'*H*B.OR; % project the block right Hamiltonian

  % project S2 into the new block
  % we don't need h2 because it's already included in H
  temp = speye(length(B.OR)/length(B.S2{1}));
  for k=1:length(B.S2)
    B.S2p{k} = B.OR'*kron(B.S2{k}, temp)*B.OR;
  end
end
