function L = superop_lindblad(A, H)
% SUPEROP_LINDBLAD  Liouvillian superoperator for a set of Lindblad operators.
%  L = superop(A [, H])
%
%  A is a cell vector of traceless, orthogonal Lindblad operators.
%  H is an optional Hamiltonian operator.
%
%  Returns the Liouvillian superoperator L corresponding to the
%  diagonal-form Lindblad equation
%
%  \dot{\rho} = inv_vec(L * vec(\rho)) =
%    = -i [H, \rho] +\sum_k (2 A_k \rho A_k^\dagger -\{A_k^\dagger A_k, \rho\})

% James D. Whitfield 2009
% Ville Bergholm 2009-2010


if (nargin == 2)
  % Hamiltonian
  iH = i*H;
else
  iH = 0;
end

L = 0;
acomm = 0;
for k=1:length(A)
  acomm = acomm + A{k}' * A{k};
  L = L +lrmul(2 * A{k}, A{k}'); 
end
L = L +lmul(-acomm -iH) +rmul(-acomm +iH);