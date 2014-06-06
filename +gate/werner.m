function rho = werner(p, d)
% WERNER  Werner states.
%  rho = werner(p [, d=2])
%
%  For every d >= 2, Werner states are a linear family of
%  bipartite d \times d  dimensional quantum states that are
%  invariant under all local unitary rotations of the form
%  U \tensor U, where U \in SU(d).
%
%  p \in [0,1] is the weight of the symmetric part of the state.
%
%  The state is entangled when p < 1/2, and pure only when d = 2
%  and p = 0, at which point it becomes the 2-qubit singlet state.

%! R.F.Werner, "Quantum states with Einstein-Podolsky-Rosen correlations admitting a hidden-variable model". PRA 40, 4277 (1989). doi:10.1103/PhysRevA.40.4277

% Ville Bergholm 2014

if nargin < 2
    d = 2;  % two qubits by default
end

S = gate.swap(d, d);
I = gate.id([d, d]);

temp = 1 -2*p;
alpha = (d*temp +1) / (temp +d);

rho = (I -alpha*S) / (d * (d -alpha));
rho = state(rho);
rho.data = full(rho.data);

% parameter value which yields I
%p_id = (d+1)/(2*d)

%d = 3; res = []; for k=1:length(p), r = werner(p(k),d); res(k,:) = [r.purity(), r.lognegativity(1)]; end; plot(p, res);
