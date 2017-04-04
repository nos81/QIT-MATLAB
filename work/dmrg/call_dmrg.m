% interface test script for finite dmrg
% Ville Berholm 2010

sz = qit.sz;

do_tdmrg = false; % tdmrg is not reliable at the moment


%====== Define the system =======

dims = 2 * ones(1, 4); % dimensions of the sites
n = length(dims) % number of sites

% Hamiltonian

%J = @(k,s) cos(0.5*pi*(k-1)/n);
%J = @(k,s) cos(0.05 + pi*(k-1)/n);
%J = @(k,s) 2*heaviside(k-n/2 +0.1)-1;
J = @(k,s) 1;

%h = @(k) 0; 
h = @(k) cos(0.05 + 2.4*pi*(k-1)/n);
% NOTE if a site happens to be uncoupled from the rest of the chain
% and has a negligible local H (within numerical precision), its
% ground state is ill-defined and thus more or less random!

Hfunc = @(ss) hamiltonian.heisenberg(dims, J, h, ss);



%====== Finite-system DMRG =======

m = 30 % truncation dimension

block = [] % DMRG block data structure
[E, block] = dmrg_finite(Hfunc, n, m, 3);
E



%====== Time-dependent DMRG =======

if (do_tdmrg)
  disp('t-DMRG...')
  % propagate the system in time using another Hamiltonian
  H2func = @(ss) hamiltonian.heisenberg(dims, J, 0, ss);
  t = 0.374
  steps = 50;
  block = tdmrg(block, H2func, t, steps, m);
end



%====== Measurements =======

% measured quantities
observables = {{sz, 1}, {sz, 2}, {sz, n-1}, {sz, n}, {sz, 1; sz, 2}}; %, {sz, 2; sz, 3}};

% measure the DMRG observables
res_dmrg = cell2mat(dmrg_measure(block, observables))


res_exact = [];
if (n <= 10)
  % if the system is small enough we can do exact diagonalization for comparison
  disp('Diagonalizing the exact Hamiltonian for comparison...')
  [H, dim] = Hfunc([1 n]);
  [v,d] = eig(full(H));
  s = state(v(:,1), dim);
  E0 = d(1,1)

  % exact t propagation
  if (do_tdmrg)
    [K, dim] = H2func([1 n]);
    s = s.prop(expm(-i*t*K));
  end

  % and measurement
  for k=1:length(observables)
    res_exact(k) = ev(s, op_list({observables{k}}, dim));
  end
  
  % subtract exact ground state energy E0 from the DMRG estimates
  E = abs(E-E0);
  
  res_exact
  res_diff = res_dmrg-res_exact
else
  % we have no exact E0 value so just fix the lowest obtained energy to 1e-20 for the semilog plot
  E = 1e-20 +E-min(min(E));
end

figure
semilogy(E.')
temp = sprintf('Energy convergence in DMRG sweeps, %d sites, m = %d', n, m);
title(temp)
xlabel('site')
ylabel('E-E_{min}')




% questions

% we do not sweep the end cells. why is this ok?
% because for these end blocks dim < m and every detail is kept?
% after the initial chain growth sweep we need not visit them again?
% the unitary O only amounts to a rotation of the observables local to these blocks?
% if so, can we ignore more sites from both ends as long as dim < m?
% - yes, it seems so

% at ends of the sweeps, do we do the same diagonalization twice??!!?
% is it enough to do it only once?
% - yes


% idea:
% save some CPU by skipping the ends in the sweeps after initial chain growth sweep.
% the end-block H's will be OK, but the end-block states will not
% (they will correspond to the growth-phase fake state). before
% measurements or time propagation, re-compute the end block states.
