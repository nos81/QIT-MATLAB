function [E0, Et] = dmrg_finite(H_func, n, m, sweeps)
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

% precompute ops (if constant, this is a waste)
for k=1:n-1
  [block{k}.S1, block{k}.S2, block{k}.h1, block{k}.h2] = H_func(k);    
end


% block k:   L(1:k)|R(k+1:n), S1|S2
% contents: HL, HR
% vector S1: acting on the rightmost site represented by the L block
% vector S2: acting on the leftmost site represented by the R block
% h1, h2 likewise

b = 1 % initial block size

% initial L block
block{b}.P = 0;
H = H_func([1 b]);
block{b} = projectL(block{b}, speye(length(H)), H);

% initial R block
block{n-b}.P = 0;
H = H_func([n-b+1 n]);
block{n-b} = projectR(block{n-b}, speye(length(H)), H);


% grow chain
%if (false)
if (true)

% TODO how about L*|R ?
    
% L*|*R scheme (heuristic)
for k=b+1:ceil(n/2)
  % build block Hamiltonians
  prev = k-1;
  [HL, Ldim] = Lstar(block{prev});
  [HR, Rdim] = starR(block{n-prev});
  
  % superblock Hamiltonian
  % NOTE left coupling
  % the ** coupling is undefined until the blocks touch so we just make it up!
  H = kron(HL, speye(length(HR))) +kron(speye(length(HL)), HR)...
      +coupling(speye(Ldim), block{k}.S1, block{k}.S2, speye(Rdim));
  H = 0.5*(H+H'); % eliminate rounding errors

  % truncation
  [OL, OR, block{k}.P, E0] = truncate(H, [length(HL), length(HR)], m);
  E0

  block{k} = projectL(block{k}, OL, HL);
  block{n-k} = projectR(block{n-k}, OR, HR);
end

% half sweep to right
sweep(k+1:n-b-1, true);

else

% L*|**** scheme
ns = 4; % number of sites in the right block

for k=b+1:n-b-1  %n-1
  % build block Hamiltonians
  % the L* block has k sites in it
  prev = k-1;
  [HL, Ldim] = Lstar(block{prev});
  
  sites = min(ns, n-k);
  [HR, temp] = H_func([k+1, k+sites]);
  Rdim = length(HR)/temp(1);
  
  % superblock Hamiltonian
  H = kron(HL, speye(length(HR))) +kron(speye(length(HL)), HR)...
      +coupling(speye(Ldim), block{k}.S1, block{k}.S2, speye(Rdim));
  H = 0.5*(H+H'); % eliminate rounding errors

  % truncation
  [OL, OR, block{k}.P, E0] = truncate(H, [length(HL), length(HR)], m);
  E0

  block{k} = projectL(block{k}, OL, HL);
  %block{n-k} = projectR(block{n-k}, OR, HR);
  % FIXME more here?
end
end



% actual sweeps
range = b+1:n-b-1;
for s=1:sweeps
  sweep(range, false);
  sweep(range, true);
end


% accurate energies for comparison
[H, dim] = H_func([1 n]);
Et = eigs(H, 1, 'SA'); % gs energy
E0-Et


function sweep(sites, sweep_right)
  if (~sweep_right)
    sites = fliplr(sites);
  end

  for k=sites
    % after this iteration, block k will be replaced

    % build block Hamiltonians for L(k-1)* and *R(k+1)
    [HL, Ldim] = Lstar(block{k-1});
    [HR, Rdim] = starR(block{k+1});
    
    % superblock Hamiltonian L(k-1)*|*R(k+1)
    H = kron(HL, speye(length(HR))) +kron(speye(length(HL)), HR)...
        +coupling(speye(Ldim), block{k}.S1, block{k}.S2, speye(Rdim));
    H = 0.5*(H+H'); % eliminate rounding errors

    % truncation
    [OL, OR, block{k}.P, E0] = truncate(H, [length(HL), length(HR)], m);
    E0

    if (sweep_right)
      block{k} = projectL(block{k}, OL, HL);
    else
      block{k} = projectR(block{k}, OR, HR);
    end
  end  
end

function [H, Ldim] = Lstar(B)
% L* block
  Ldim = length(B.HL);
  sdim = length(B.h2);
  H = kron(B.HL, speye(sdim)) +coupling(1, B.S1p, B.S2, 1) +kron(speye(Ldim), B.h2);
end

function [H, Rdim] = starR(B)
% *R block
  sdim = length(B.h1);
  Rdim = length(B.HR);
  H = kron(speye(sdim), B.HR) +coupling(1, B.S1, B.S2p, 1) +kron(B.h1, speye(Rdim));
end

end



function [OL, OR, P, E0] = truncate(HS, dim, m)
  % find ground state
  [v,E0] = eigs(HS, 1, 'SA'); % 'SR'

  % Schmidt decompose into two halves
  s = state(v, dim);
  [d, u, v] = schmidt(s, [1]);

  keep = min(m, length(d)); % how many states to keep in A?
    
  P = 1 - sum(d(1:keep).^2); % truncation error
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
