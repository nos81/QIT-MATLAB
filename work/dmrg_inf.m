function [E, P, Et] = dmrg_inf(N, m, sdim)
% infinite spin chain, DMRG method test

%@Article{PhysRevLett.69.2863,
%  title = {Density matrix formulation for quantum renormalization groups},
%  author = {White, Steven R.},
%  journal = {Phys. Rev. Lett.},
%  volume = {69},
%  pages = {2863--2866},
%  year = {1992},
%  doi = {10.1103/PhysRevLett.69.2863}}

% Ville Bergholm 2010


if (nargin < 3)
  sdim = 2; % spin dimension
  if (nargin < 2)
    m = 5; % states to keep
  end
end

J = [1 1 1];
h = @(k) 0;

b = 1 % initial block size

E = zeros(1,b);
P = zeros(1,b);

% initial block Hamiltonian
[HA, dim] = hamiltonian.heisenberg(sdim*ones(1,b), J, h);

I_S = speye(sdim);

% initial ops
[S1, S2, ~, h2] = hamiltonian.heisenberg([sdim sdim], J, h, b);
SA = project_S(S1, speye(length(HA))); % kinda wasteful but

if (true)
  % A**A scheme, infinite chain
  for len=b+1:N
    % after this iteration, block A will represent len sites

    dimA = length(HA); % m, except on the first rounds...
    I_A = speye(dimA);
    
    % build block Hamiltonian for A*
    HB = kron(HA, I_S) +coupling(SA, 1, S2) +kron(I_A, h2);

    % Hamiltonian operators for sites {len, len+1} (**)
    [S1, S2, ~, h2] = hamiltonian.heisenberg([sdim sdim], J, h, len);

    % superblock Hamiltonian
    % clever? A*|A*, couple the * sites using non-nn topology, no reflection
    I_B = speye(dimA*sdim);
    HS = kron(HB, I_B) +kron(I_B, HB) +coupling(S1, I_A ,S2);
    HS = 0.5*(HS+HS'); % eliminate rounding errors
    
    % find ground state
    [v,d] = eigs(HS, 1, 'SA'); % 'SR'
    E(len) = d/(2*len); % energy/site

    % Schmidt decompose into two halves
    s = state(v, [dimA sdim dimA sdim]);
    [d, u, v] = schmidt(s, [1 2]);

    keep = min(m, size(u, 1)); % how many states to keep in A?
    
    P(len) = 1 - sum(d(1:keep).^2); % truncation error
    O{len} = u(:, 1:keep); % keep lowest eigenstates
    
    HA = O{len}'*HB*O{len}; % project the block Hamiltonian
    
    % project S1 to get SA for the next round
    % we don't need h1 because it's already included in HA
    SA = project_S(S1, O{len});
  end
  
else
  % A*.(times k) scheme, infinite chain
  k = 4;
  
  % k new sites H
  [H2, dim] = hamiltonian.heisenberg(sdim*ones(1,k), J, h);
  dim2 = length(H2);
  I_2 = speye(dim2);
  
  for len=b+1:N
    dimA = length(HA); % m, except on the first rounds...
    I_A = speye(dimA);

    % new block H with just one extra site
    HB = kron(HA, I_S) +coupling(SA, 1, S2) +kron(I_A, h2);
    I_B = speye(dimA*sdim);
    
    % Hamiltonian operators for sites {len, len+1} (**)
    [S1, S2, ~, h2] = hamiltonian.heisenberg([sdim sdim], J, h, len);
    
    % superblock: add the k new sites
    HS = kron(HB, I_2) +mkron(I_A, coupling(S1, 1, S2), I_rest) +kron(I_B, H2);
    HS = 0.5*(HS+HS'); % eliminate rounding errors

    % find ground state
    [v,d] = eigs(HS, 1, 'SA'); % 'SR'
    E(len) = d/(len+k); % energy/site

    % Schmidt decompose into two halves
    s = state(v, [dimA, sdim*ones(1, 1+k)]);
    [d, u, v] = schmidt(s, [1 2]);

    keep = min(m, size(u, 1)); % how many states to keep in A?
    
    P(len) = 1 - sum(d(1:keep).^2); % truncation error
    O{len} = u(:, 1:keep); % keep lowest eigenstates
    
    HA = O{len}'*HB*O{len}; % project the block Hamiltonian
    
    % project S1 to get SA for the next round
    % we don't need h1 because it's already included in HA
    SA = project_S(S1, O{len});
  end
end



% accurate energies for comparison
sss = 1:10;
for k=1:length(sss)
  [H, dim] = hamiltonian.heisenberg(sdim*ones(1,sss(k)), J, h);
  Et(k) = eigs(H, 1, 'SA')/sss(k); % gs energy per spin
end

E0 = 0.25 -log(2); % exact answer (Bethe Ansatz) for infinite spin-1/2 chain, Heisenberg, antiferromagnetic
E = E-E0;
Et = Et-E0;

figure
hold on
semilogy(Et, 'k-')
%plot(Et, 'k-')
semilogy(E, 'r-')
%plot(E, 'r-')
xlabel('sites')
ylabel('(E-E_0)/site')

%max(nh)

HA
end



function SA = project_S(S, O)
%
  temp = speye(size(O, 1)/size(S{1}, 1));
  for k=1:length(S)
    SA{k} = O'*kron(temp, S{k})*O;
  end
end


function C = coupling(S1, I_A, S2)
% C = \sum_k I_A * S1_k * I_A * S2_k
  C = sparse(0);
  for k=1:length(S1)
    C = C +mkron(S1{k}, I_A, S2{k});
  end
  C = kron(I_A, C);
end
