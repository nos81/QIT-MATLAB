function [H_low, H_eff, Q] = gadget(type, A, B, C, H_rest)
% GADGET  Perturbative gadgetry.
%  [H_low, H_eff, Q] = gadget('kkr_3to2', A, B, C, H0) % Kempe-Kitaev-Regev 3-to-2 gadget
%  [H_low, H_eff, Q] = gadget('ot_3to2',  A, B, C, H0) % Oliveira-Terhal 3-to-2 gadget
%  [H_low, H_eff, Q] = gadget('ot_sd',    A, B, H0)    % Oliveira-Terhal subdivision gadget
%
%  Builds several different pertubative Hamiltonian gadgets for approximating k-local
%  Hamiltonians with Hamiltonians of a lower degree of locality.
%
%  The output variables are H_low, the restriction of the gadget
%  Hamiltonian to the low eigenspace, H_eff, the effective
%  Hamiltonian, and Q = H+V, the full gadget Hamiltonian.
%
%  The effective Hamiltonians for different gadgets are
%  kkr_3to2: H_eff = H0 -6 * A \otimes B \otimes C
%  ot_3to2:  H_eff = H0 +A \otimes B \otimes C
%  ot_sd:    H_eff = H0 +A \otimes B

%! J. Kempe, A. Kitaev and O. Regev, "The Complexity of the Local Hamiltonian Problem", SIAM Journal of Computing 35, 1070 (2006). arXiv.org:quant-ph/0406180
%! R. Oliveira and B.M. Terhal, "The complexity of quantum spin systems on a two-dimensional square lattice", Quant. Inf. Comp. 8, 0900 (2008). arXiv.org:quant-ph/0504050
% Ville Bergholm 2009


if (nargin < 1)
  error('Specify gadget type.');
end

switch lower(type)
  case 'kkr_3to2'
  case 'ot_3to2'
    if (nargin < 5)
      H_rest = 0;
      if (nargin < 4)
        error('Need A, B, and C.');
      end
    end

  case 'ot_sd' 
    if (nargin < 4)
      H_rest = 0;
      if (nargin < 3)
        error('Need A and B.');
      end
    else
      H_rest = C;
    end
    C = 0; % just for convenience
end


epsilon = 0.001;

global qit;

% mediator qubit ops
W_I = qit.I;
W_0 = qit.p0;
W_1 = qit.p1;
W_X = qit.sx;
W_Z = qit.sz;


r = max([norm(A), norm(B), norm(C)])

dim = [size(A, 1), size(B, 1), size(C, 1)]
d_ops = prod(dim);

A = mkron(A, eye(prod(dim(2:3))));
B = mkron(eye(dim(1)), B, eye(dim(3)));
C = mkron(eye(prod(dim(1:2))), C);

switch lower(type)
  case 'kkr_3to2'
    % kempe-kitaev-regev 3-to-2 local gadget
    % A,B,C must be positive semidefinite (nonneg. eigenvalues)
    % uses three mediator qubits
    Delta = epsilon^(-3); % TODO should be a function of r too

    H_target = H_rest -6*A*B*C;

    H = kron(-Delta/4*(mkron(W_Z, W_Z, W_I) +mkron(W_Z, W_I, W_Z) +mkron(W_I, W_Z, W_Z) -3*eye(8)), eye(d_ops));
    X = H_rest +Delta^(1/3)*(A^2 +B^2 +C^2);
    V = kron(eye(8), X) -Delta^(2/3)*(mkron(W_X, W_I, W_I, A) +mkron(W_I, W_X, W_I, B) +mkron(W_I, W_I, W_X, C));

  case 'ot_sd'
    % subdivision gadget
    H_restp = H_rest +A^2/2 +B^2/2;
    H_target = H_restp -(-A+B)^2/2; % == H_rest +A*B;

    assert(max(norm(H_restp), r) >= 1)
    % Delta must be >= 1 ???

    Delta = (norm(H_restp) + sqrt(2)*r)^6/epsilon^2

    H = kron(Delta*W_1, eye(d_ops));
    V = kron(W_I, H_restp) +kron(sqrt(Delta/2)*W_X, -A+B);

  case 'ot_3to2'
    % 3-local to 2-local gadget
    Delta = epsilon^(-3); % TODO should be a function of r

    H_target = H_rest +A*B*C;

    H = kron(Delta*W_1, eye(d_ops));

    V_extra = Delta^(1/3)*(-A+B)^2/2 +(A^2 + B^2)*C/2;
    V = kron(W_I, H_rest+V_extra) +kron(-Delta^(2/3)*W_1, C) +kron(Delta^(2/3)*W_X/sqrt(2), -A+B);
end

Q = H+V;

%H_eff = kron(W_0, H_target);
H_eff = H_target;

% construct the projector P_-
[v,d] = eig(Q);
d = diag(d);
ind = find(d < Delta/2);
d(ind)
P = v(:,ind)*v(:,ind)';
%P(find(abs(P) < 1e-5)) = 0;

H_low = P*Q*P;
H_low = (H_low+H_low')/2; % eliminate numerical errors, it should be hermitian
% H_low(1:16,1:16) ~== H_eff, why???
