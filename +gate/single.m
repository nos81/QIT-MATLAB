function U = single(L, targ, dim)
% SINGLE  Single-qudit operator.
%
%  U = single(L, targ, dim)
%
%  Returns the operator U corresponding to the local operator L applied
%  to subsystem targ (and identity applied to the remaining subsystems).
%
%  dim is either the dimension vector for U or an integer scalar denoting
%  the number of subsystems in an all-qubit system.

% James Whitfield 2010
% Ville Bergholm 2010


if (isscalar(dim))
  dim = 2*ones(1,dim); % assume qubits
end

U = op_list({{L, targ}}, dim);
