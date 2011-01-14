function T = tangle(s, sys)
% TANGLE  Tangle of the state.
%  T = tangle(s, sys)
%
%  Returns the tangle of the state s wrt. the partitioning
%  given by the listing of subsystems in the vector sys.

%! V. Coffman, J. Kundu, W.K. Wootters, "Distributed entanglement", PRA 61, 052306 (2000). doi: 10.1103/PhysRevA.61.052306
% Ville Bergholm 2010

% 3-tangle only, qubits

% 2-tangle is squared concurrence
% T3 = T2(A:BC) -T2(A:B) -T2(A:C)

%A = ptrace(s, [2 3]);
%[x,d] = eig(A.data)
%diag(d)
%[lambda, u, v] = schmidt(s, 1)
%lambda.^2
%abs(v(:,2)'*x(:,3))

AB = ptrace(s, [3]);
AC = ptrace(s, [2]);

T = concurrence(s, [1])^2 -concurrence(AB)^2 -concurrence(AC)^2;
