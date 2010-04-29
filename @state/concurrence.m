function C = concurrence(s, sys)
% CONCURRENCE  Concurrence of the state.
%  C = concurrence(s, sys)
%
%  Returns the concurrence of the state s wrt. the partitioning
%  given by the listing of subsystems in the vector sys.

%! W.K. Wootters, "Entanglement of Formation of an Arbitrary State of Two Qubits", PRL 80, 2245 (1998).
%! R. Horodecki, P. Horodecki, M. Horodecki, K. Horodecki, "Quantum entanglement", arXiv:quant-ph/0702225 (2007).
% Ville Bergholm 2006-2008 




% FIXME
% bipartite:
% entropy of entanglement E: -tr(rho_a * log2(rho_a))
% e of formation: minimum of average E over pure decompositions
% distillable e, e cost, relative entropy of e.

sy = [0 -i; i 0];
X = kron(sy,sy);

p = s.data;
if (size(p, 2) == 1)
  % ket
  p = p/norm(p); % normalisoidaan tila
  C = abs(p'*X*conj(p));

  rho_A = ptrace(s, sys);
  temp = rho_A.data;
  xxx = sqrt(2*(1-real(trace(temp*temp))))

  % konkurrenssi, p tilavektori normaalikannassa, a sama magic basessa

  % phi+-, psi+- eli tripletti-singletti
  bell = [1 i 0 0; 0 0 i 1; 0 0 i -1; 1 -i 0 0]/sqrt(2);
  a = bell'*p;
  yyy = abs(sum(a.^2))

else
  % state operator
  p = p/trace(p); % normalisoidaan

  if (abs(trace(p*p)-1) > 1e-6)
    disp('not a pure state!');
    l = sqrt(sort(real(eig(p*X*conj(p)*X)), 'descend'))
    C = max(0, l(1)-l(2)-l(3)-l(4));
  else
    C = sqrt(real(trace(p*X*conj(p)*X))); % same formula as for state vecs
  end
end
