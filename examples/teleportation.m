function [reg_B, payload] = teleportation()
% TELEPORTATION  Quantum teleportation demo.
%
%  [reg_B, payload] = teleportation()
%
%  Simulate the teleportation of a qubit from Alice to Bob.

% Ville Bergholm 2009


fprintf('\n\n=== Quantum teleportation ===\n\n')

global qit;

H    = gate.walsh(1); % Hadamard gate
cnot = gate.controlled(qit.sx, 1);
I    = qit.I;

disp('Alice and Bob start with a shared EPR pair.')
epr = state('bell1')


disp('Alice wants to transmit this payload qubit to Bob:')
payload = state('0');
payload = u_propagate(payload, rand_SU(2));
% choose a nice global phase
payload = fix_phase(payload)


disp('The total |payload> \otimes |epr> register looks like')
reg = tensor(payload, epr)


disp('Now Alice entangles the payload with her half of the EPR pair')
reg = u_propagate(reg, kron(H, eye(4)) * kron(cnot, I))

[p, b(1), reg] = measure(reg, 1);
[p, b(2), reg] = measure(reg, 2);
%[p, b, reg] = measure(reg, [1 2]);
b = b-1;
fprintf('and measures her qubits, getting the result [%d, %d].\n', b)
disp('She then transmits the two bits to Bob using a classical channel. The shared state is now')
reg

disp('Since Alice''s measurement has unentangled the state,')
disp('Bob can ignore her qubits. His qubit now looks like')
reg_B = to_ket(ptrace(reg, [1 2])) % pure state

disp('Using the two classical bits of data Alice sent him,')
disp('Bob performs a local transformation on his half of the EPR pair.')
reg_B = fix_phase(u_propagate(reg_B, qit.sz^(b(1)) * qit.sx^(b(2))))


ov = overlap(payload, reg_B);
fprintf('The overlap between the resulting state and the original payload state is |<payload|B>| = %f\n', ov)
if (norm(ov-1) > qit.tol)
  error('Should not happen.')
else
  disp('The payload state was succesfully teleported from Alice to Bob.')
end
