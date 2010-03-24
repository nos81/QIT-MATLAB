function [E] = dmrg_finite(H_func, n, m, sweeps)
% DMRG_FINITE  DMRG method for solving finite 1d-nn quantum systems.
%  [E0] = dmrg_finite(H_func, n, m, sweeps)
%
%  Returns the approximate ground state energy for a
%  one-dimensional chain of n sites with nearest-neighbor
%  couplings, with a Hamiltonian defined by H_func, using a
%  finite-system density matrix renormalization group method.
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
block{b} = projectL(block{b}, speye(length(H)), H);

% initial R block
block{n-b}.P = 0;
H = H_func([n-b+1 n]);
block{n-b} = projectR(block{n-b}, speye(length(H)), H);


% grow the chain
E = zeros(1,n-1);

if (1)
  % L*|*R scheme (heuristic)
  for k=b+1:ceil(n/2)-1
    % build block Hamiltonians
    prev = k-1;
    [HL, Ldim] = Lstar(block{prev});
    [HR, Rdim] = starR(block{n-prev});

    % NOTE the ** coupling is undefined until the blocks touch so we just make it up!
    % The superblock Hamiltonian is formed using the left coupling only.
    [block{k}, OL, OR, E(k)] = combine_and_truncate(block{k}, HL, HR, [Ldim Rdim], m);

    % discard useless states (in the sweep phase the dims are different)
    block{k} = rmfield(block{k}, 'state');
    
    block{k}   = projectL(block{k},   OL, HL);
    block{n-k} = projectR(block{n-k}, OR, HR);
  end

  % half sweep to right
  E(end+1,:) = sweep(k+1:n-b-1, true);
  
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
  
    [block{k}, OL, OR, E(k)] = combine_and_truncate(block{k}, HL, HR, [Ldim Rdim], m);

    % discard useless states (in the sweep phase the dims are different)
    block{k} = rmfield(block{k}, 'state');

    block{k} = projectL(block{k}, OL, HL);
    % projectR is unnecessary here
  end
end


% actual sweeps
range = b+1:n-b-1; % FIXME one more site?
for s=1:sweeps
  E(end+1,:) = sweep(range, false);
  E(end+1,:) = sweep(range, true);
end

% measurement sweep
%result = {};
%sweep(range, true, M);



function res = sweep(sites, sweep_right, M)
% L*|*R sweeping
% TODO how about L*|R ?

  if (~sweep_right)
    sites = fliplr(sites);
  end

  res = zeros(1,n-1);
  
  for k=sites
    % after this iteration, block k will be replaced

    % build block Hamiltonians for L(k-1)* and *R(k+1)
    [HL, Ldim] = Lstar(block{k-1});
    [HR, Rdim] = starR(block{k+1});
    
    [block{k}, OL, OR, res(k)] = combine_and_truncate(block{k}, HL, HR, [Ldim Rdim], m);
    
    if (sweep_right)
      block{k} = projectL(block{k}, OL, HL);
    else
      block{k} = projectR(block{k}, OR, HR);
    end
  end  
end

end


function [H, dim] = Lstar(B)
% L* block
  dim = [length(B.HL), length(B.h2)];
  H = kron(B.HL, speye(dim(2))) +coupling(1, B.S1p, B.S2, 1) +kron(speye(dim(1)), B.h2);
end

function [H, dim] = starR(B)
% *R block
  dim = [length(B.h1), length(B.HR)];
  H = kron(speye(dim(1)), B.HR) +coupling(1, B.S1, B.S2p, 1) +kron(B.h1, speye(dim(2)));
end

function [B, OL, OR, E0] = combine_and_truncate(B, HL, HR, dim, m)

  % superblock Hamiltonian
  H = kron(HL, speye(length(HR))) +kron(speye(length(HL)), HR)...
      +coupling(speye(dim(1)), B.S1, B.S2, speye(dim(4)));

  H = 0.5*(H+H'); % eliminate rounding errors

  % find target state
  opts = struct();
  if (isfield(B, 'state'))
    opts.v0 = B.state.data; % use last round's result as a guess
  end
  [v,E0] = eigs(H, 1, 'SA', opts); % ground state
  B.state = state(v, dim); % store it
  
  % Schmidt decompose into two halves
  s = state(v, [length(HL), length(HR)]);
  [d, u, v] = schmidt(s, [1]);

  keep = min(m, length(d)); % how many states to keep in A?
    
  B.P = 1 - sum(d(1:keep).^2); % truncation error
  % keep most influential states
  OL = u(:, 1:keep);
  OR = v(:, 1:keep);
end


function B = projectL(B, O, H)
  B.HL = O'*H*O; % project the block left Hamiltonian
  
  % project S1 into the new block
  % we don't need h1 because it's already included in H
  temp = speye(length(O)/length(B.S1{1}));
  for k=1:length(B.S1)
    B.S1p{k} = O'*kron(temp, B.S1{k})*O;
  end
end

function B = projectR(B, O, H)
  B.HR = O'*H*O; % project the block right Hamiltonian

  % project S2 into the new block
  % we don't need h2 because it's already included in H
  temp = speye(length(O)/length(B.S2{1}));
  for k=1:length(B.S2)
    B.S2p{k} = O'*kron(B.S2{k}, temp)*O;
  end
end


function C = coupling(IA, S1, S2, IB)
% C = \sum_k IA * S1_k * IB * S2_k
  C = sparse(0);
  for k=1:length(S1)
    C = C +kron(S1{k}, S2{k});
  end
  C = mkron(IA, C, IB);
end
