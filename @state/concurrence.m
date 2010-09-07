function C = concurrence(s, sys)
% CONCURRENCE  Concurrence of the state.
%  C = concurrence(s, sys)
%
%  Returns the concurrence of the state s wrt. the partitioning
%  given by the listing of subsystems in the vector sys.

%! W.K. Wootters, "Entanglement of Formation of an Arbitrary State of Two Qubits", PRL 80, 2245 (1998).
%! R. Horodecki, P. Horodecki, M. Horodecki, K. Horodecki, "Quantum entanglement", arXiv:quant-ph/0702225 (2007).
% Ville Bergholm 2006-2010


global qit;

if (abs(trace(s)-1) > qit.tol)
  disp('Warning: State not properly normalized.')
end

if (nargin == 2)
  % concurrence between a qubit and a larger system

  if (length(sys) == 1 && s.dim(sys) == 2)
    if (abs(purity(s)-1) < qit.tol)
      % pure state

      n = length(s.dim);
      rho_A = ptrace(s, setdiff(1:n, sys)); % trace over everything but sys
      C = 2 * sqrt(real(det(rho_A.data))); % = sqrt(2*(1-real(trace(temp*temp))))
      return
    end
    error('Not a pure state.')

  end
  error('Concurrence only defined between a qubit and another system.')
end


if (length(s.dim) ~= 2 || any(s.dim ~= [2 2]))
  % not a two-qubit state
  error('Not a two-qubit state.')
end


X = kron(qit.sy, qit.sy);

p = s.data;
if (size(p, 2) == 1)
  % ket
  C = abs(p'*X*conj(p));

  % find the coefficients a of the state ket in the magic base
  % phi+-, psi+-,  = triplet,singlet
  %bell = [1 i 0 0; 0 0 i 1; 0 0 i -1; 1 -i 0 0]/sqrt(2);
  %a = bell'*p;
  %C = abs(sum(a.^2))

else
  % state operator

  if (abs(purity(s)-1) > qit.tol)
    l = real(sqrt(sort(real(eig(p * X * conj(p) * X)), 'descend')));
    C = max(0, l(1) -l(2) -l(3) -l(4));
  else
    C = sqrt(real(trace(p*X*conj(p)*X))); % same formula as for state vecs
  end
end
