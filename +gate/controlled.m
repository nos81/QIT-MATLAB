function out = controlled(U, ctrl)
% CONTROLLED  Controlled-unitary gate.
%  V = controlled(U, ctrl)
%
%  Returns the (t+n)-qubit controlled-U gate, where t == length(ctrl) and n == size(U,1).
%  ctrl is a vector defining the control nodes. It has one entry per control qubit.
%  0 denotes control on zero, 1 denotes control on one, and -1 denotes no control.
%
%  controlled(NOT, 1) gives the standard CNOT gate. controlled(NOT, [1 1]) gives the Toffoli gate.

% Ville Bergholm 2009
% TODO enable arbitrary permutation of control and target qubits?


% no control, control on zero, control on one
c = {[1 1], [1 0], [0 1]} ;

t = length(ctrl);
pass = 1;
for k=1:t
    pass = kron(pass, c{ctrl(k)+2});
end

fail = ones(1, 2^t) - pass;
n = size(U,1);
out = diag(kron(fail, ones(1,n))) + kron(diag(pass), U);
