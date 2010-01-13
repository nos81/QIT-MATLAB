function superdense_coding()
% SUPERDENSE_CODING  Superdense coding demo.
%
%  superdense_coding()
%
%  Simulate Alice sending two bits of information to Bob using superdense coding.

% Ville Bergholm 2010


fprintf('\n\n=== Superdense coding ===\n\n')

global qit;

H    = qit.H; % Hadamard gate
cnot = gate.controlled(qit.sx, 1);
I    = qit.I;

disp('Alice and Bob start with a shared EPR pair.')
reg = state('bell1')

% two random bits
a = round(rand(1, 2));
fprintf('Alice wishes to send two bits of information to Bob: a = [%d, %d].\n', a)

disp('Alice encodes the bits to her half of the EPR pair using local transformations,')
reg = u_propagate(reg, kron(qit.sz^(a(1)) * qit.sx^(a(2)), I))

disp('and sends it to Bob. He then disentangles the pair using a CNOT and a Hadamard,')
reg = u_propagate(reg, kron(H, I) * cnot)

[p, b(1), reg] = measure(reg, 1);
[p, b(2), reg] = measure(reg, 2);
%[p, b, reg] = measure(reg, [1 2]);
b = b-1;
fprintf('and measures both qubits in the computational basis, obtaining the result  b = [%d, %d].\n', b)

if (a == b)
  disp('The bits were transmitted succesfully.')
else
  error('Should not happen.')
end
