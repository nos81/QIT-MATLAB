function [reg_B, payload] = teleportation()
% TELEPORTATION  Quantum teleportation.
%
%  [reg_B, payload] = teleportation()
%
%  Simulate the teleportation of a qubit from Alice to Bob.

% Ville Bergholm 2009


fprintf('\n\n=== Quantum teleportation ===\n')

global qit;

H    = gate.walsh(1); % Hadamard gate
cnot = gate.controlled(qit.sx, 1);
I    = qit.I;

disp('Alice and Bob start with a shared EPR pair.')
epr = state('bell1')

disp('Alice wants to transmit this payload qubit to Bob:')
payload = state('0');
payload.data = rand_SU(2)*payload.data;
% choose a nice global phase
temp = payload.data;
phase = temp(1)/abs(temp(1));
payload.data = payload.data/phase;
payload

disp('The total |payload> \otimes |epr> register looks like')
reg = tensor(payload, epr)


disp('Now Alice entangles the payload with her half of the EPR pair')
reg.data = kron(H, eye(4)) * kron(cnot, I) * reg.data

disp('and then measures her qubits, getting the results')
%[p, q, reg] = measure(reg, [1 2]);
[p, q1, reg] = measure(reg, {mkron(qit.p0, I, I), mkron(qit.p1, I, I)});
[p, q2, reg] = measure(reg, {mkron(I, qit.p0, I), mkron(I, qit.p1, I)});
q1 = q1-1
q2 = q2-1
disp('and transmits them to Bob using a classical channel. The shared state is now')
reg

disp('Since Alice''s measurement has unentangled the state,')
disp('Bob can ignore her qubits. His qubit now looks like')
reg_B = to_ket(ptrace(reg, [1 2])) % pure state

disp('Using the two classical bits of data Alice sent him,')
disp('Bob performs a local transformation on his half of the EPR pair.')
%reg.data = kron(eye(4), qit.sz^(q1) * qit.sx^(q2)) * reg.data
reg_B.data = qit.sz^(q1) * qit.sx^(q2) * reg_B.data


ov = overlap(payload, reg_B);
fprintf('The overlap between the resulting state and the original payload state is |<payload|B>|^2 = %f\n', ov)
if (norm(ov-1) > qit.tol)
  error('Should not happen.')
else
  disp('The payload state was succesfully teleported from Alice to Bob.')
end
