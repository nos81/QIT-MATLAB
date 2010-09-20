function out = controlled(U, ctrl, dim)
% CONTROLLED  Controlled gate.
%  V = controlled(U, ctrl [, dim])
%
%  Returns the (t+1)-qudit controlled-U gate, where t == length(ctrl).
%  ctrl is an integer vector defining the control nodes. It has one entry k per
%  control qudit, denoting the required computational basis state |k>
%  for that particular qudit. Value k == -1 denotes no control.
%
%  dim is the dimensions vector for the control qudits. If not given, all controls
%  are assumed to be qubits.
%
%  Examples:
%   controlled(NOT, [1]) gives the standard CNOT gate.
%   controlled(NOT, [1 1]) gives the Toffoli gate.

% Ville Bergholm 2009-2010


t = length(ctrl);

if (nargin < 3)
  dim = 2*ones(1, t); % qubits by default
end

if (t ~= length(dim))
  error('ctrl and dim vectors have unequal lengths.')
end

if (any(ctrl >= dim))
  error('Control on non-existant state.')
end


pass = sparse(1);
for k=1:t
  % just the diagonal
  if (ctrl(k) >= 0)
    temp = sparse(dim(k), 1);
    temp(ctrl(k)+1) = 1;
    pass = kron(pass, temp); % control on k
  else
    pass = kron(pass, ones(dim(k), 1)); % no control on this qudit
  end
end

fail = sparse(1) - pass;
T = prod(dim);
n = size(U,1);
out = kron(spdiags(fail, 0, T, T), speye(n)) + kron(spdiags(pass, 0, T, T), U);
