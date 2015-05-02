function rho = werner(p, d)
% WERNER  Werner states.
%  rho = werner(p [, d=2])
%
%  For every d >= 2, Werner states are a linear family of
%  bipartite d \times d  dimensional quantum states that are
%  invariant under all local unitary rotations of the form
%  U \tensor U, where U \in SU(d).
%
%  Every Werner state is a linear combination of the identity and
%  SWAP operators. Alternatively, they can be understood as convex
%  combinations of the appropriately scaled symmetric and
%  antisymmetric projectors P_sym = (\I+SWAP)/2 and P_asym = (\I-SWAP)/2:
%
%    \rho = p \frac{2 P_sym}{d(d+1)} +(1-p) \frac{2 P_asym}{d(d-1)}
%
%  p \in [0,1] is the weight of the symmetric part of the state.
%  The state is entangled iff p < 1/2, and pure only when d = 2
%  and p = 0, at which point it becomes the 2-qubit singlet state.
%
%  For every d, the Werner family of states includes the fully
%  depolarized state \frac{\I}{d^2}, obtained with p = \frac{d+1}{2d}.
%
%  The Werner states are partial-transpose dual to isotropic states.

%! R.F.Werner, "Quantum states with Einstein-Podolsky-Rosen correlations admitting a hidden-variable model". PRA 40, 4277 (1989). doi:10.1103/PhysRevA.40.4277

% Ville Bergholm 2014

if nargin < 2
    d = 2;  % two qubits by default
end

dim = [d, d];
S = gate.swap(d, d);
I = gate.id(dim);

%temp = 1 -2*p;
%alpha = (d*temp +1) / (temp +d);
%rho = (I -alpha*S) / (d * (d -alpha));

rho = p * (I+S)/(d*(d+1)) +(1-p) * (I-S)/(d*(d-1));
rho.data = full(rho.data);
