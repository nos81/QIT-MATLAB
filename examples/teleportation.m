function [reg_B, payload] = teleportation(d)
% TELEPORTATION  Quantum teleportation demo.
%
%  [reg_B, payload] = teleportation(d)
%
%  Simulate the teleportation of a d-dimensional qudit from Alice to Bob.

% Ville Bergholm 2009-2010


fprintf('\n\n=== Quantum teleportation ===\n\n')

global qit;

H   = gate.qft([d 1]);    % qft (generalized Hadamard) gate
add = gate.mod_add(d, d); % modular adder (generalized CNOT) gate
I   = speye(d);

% EPR preparation circuit
U = add * kron(H, I);

disp('Alice and Bob start with a shared EPR pair.')
epr = u_propagate(state('00', [d d]), U)


disp('Alice wants to transmit this payload to Bob:')
payload = state('0', d);
payload = u_propagate(payload, rand_SU(d));
% choose a nice global phase
payload = fix_phase(payload)


disp('The total |payload> \otimes |epr> register looks like')
reg = tensor(payload, epr)


disp('Now Alice entangles the payload with her half of the EPR pair')
reg = u_propagate(reg, kron(U', I))

[p, b(1), reg] = measure(reg, 1);
[p, b(2), reg] = measure(reg, 2);
%[p, b, reg] = measure(reg, [1 2]);
b = b-1;
fprintf('and measures her qudits, getting the result [%d, %d].\n', b)
disp('She then transmits the two d-its to Bob using a classical channel. The shared state is now')
reg

disp('Since Alice''s measurement has unentangled the state,')
disp('Bob can ignore her qudits. His qudit now looks like')
reg_B = to_ket(ptrace(reg, [1 2])) % pure state

disp('Using the two classical d-its of data Alice sent him,')
disp('Bob performs a local transformation on his half of the EPR pair.')
Z = diag(sqrt(d) * H(:,b(1)+1));
X = gate.mod_inc(-b(2), d);
reg_B = fix_phase(u_propagate(reg_B, Z*X))


ov = fidelity(payload, reg_B);
fprintf('The overlap between the resulting state and the original payload state is |<payload|B>| = %f\n', ov)
if (norm(ov-1) > qit.tol)
  error('Should not happen.')
else
  disp('The payload state was succesfully teleported from Alice to Bob.')
end
