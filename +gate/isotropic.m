function rho = isotropic(p, d)
% ISOTROPIC  Isotropic states.
%  rho = isotropic(p [, d=2])
%
%  For every d >= 2, isotropic states are a linear family of
%  bipartite d \times d  dimensional quantum states that are
%  invariant under all local unitary rotations of the form
%  U \tensor U^*, where U \in SU(d).
%
%  Every isotropic state is a linear combination of the identity operator
%  and the projector to the maximally entangled state |\cup> = \sum_{k=0}^{d-1} |k, k>.
%
%    \rho = p \frac{1}{d}|\cup><\cup| +(1-p) (\I-\frac{1}{d}|\cup><\cup|)/(d^2-1)
%
%  p \in [0,1] is the weight of the cup state projector in the mixture.
%  The state is entangled iff p > 1/d, and pure and fully entangled iff p = 1.
%
%  For every d, the isotropic family of states includes the fully
%  depolarized state \frac{\I}{d^2}, obtained with p = \frac{1}{d^2}.
%
%  Isotropic states are partial-transpose dual to Werner states.

% Ville Bergholm 2014

if nargin < 2
    d = 2;  % two qubits by default
end

cup = gate.copydot(0, 2, d);
cup_proj = cup * cup' / d;
I = gate.id([d, d]);

rho = p * cup_proj +(1-p) * (I -cup_proj) / (d^2 -1);
rho.data = full(rho.data);
