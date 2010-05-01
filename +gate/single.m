function U = single(g, targ, dim)
% SINGLE  gives the unitary for a single qubit gate
%
%  U = 
%
%  Returns unitary gate for a single qubit applied to target
%  qubit and identity applied to the remaining qubits. 
%
%  Counting of the qubits begins at 1. 

% James Whitfield 2010
% Ville Bergholm 2010


if (isscalar(dim))
  dim = 2*ones(1,dim); % assume qubits
end

U = op_list({{g, targ}}, dim);
