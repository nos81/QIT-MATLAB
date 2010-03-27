% interface test script for finite dmrg

blo = []
dims = 2 * ones(1, 4);

n = length(dims)

% unsymmetric chain!
%J = @(k,s) cos(0.5*pi*(k-1)/n);
J = @(k,s) cos(pi*(k-1)/n);
%J = @(k,s) 2*heaviside(k-n/2 +0.1)-1;
%J = @(k,s) -1;

%h = @(k) 0; 
h = @(k) cos(2.1*pi*(k-1)/n);
% BUG dmrg_finite randomly messes up sz(n)
% (also the state(n-1)) if n == 4


Hfunc = @(ss) hamiltonian.heisenberg(dims, J, h, ss);

m = 10
[E, blo] = dmrg_finite(Hfunc, n, m, 3);


oplist = {{qit.sz, 1}, {qit.sz, n}, {qit.sz, 1; qit.sz, 2}}; %, {qit.sz, 2; qit.sz, 3}};

% another H for propagation
H2func = @(ss) hamiltonian.heisenberg(dims, J, 0, ss);
t = 2.174
blo = tdmrg(blo, H2func, t, m)

if (n <= 6)
  % accurate energies for comparison
  [H, dim] = Hfunc([1 n]);

  [v,d] = eig(full(H));
  s = state(v(:,1), dim);
  Et = d(1,1);

  [K, dim] = H2func([1 n]);
  s = u_propagate(s, expm(-i*t*K));
  
  for k=1:length(oplist)
    sss(k) = ev(s, op_list({oplist{k}}, dim));
  end
  E = abs(E-Et);
else
  E = 1e-20+E-min(min(E));
end

%figure
%semilogy(E.')

res = dmrg_measure(blo, oplist)
sss



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
