function superdense_coding(d)
% SUPERDENSE_CODING  Superdense coding demo.
%
%  superdense_coding(d)
%
%  Simulate Alice sending two d-its of information to Bob using 
%  a shared EPR qudit pair.

% Ville Bergholm 2010


fprintf('\n\n=== Superdense coding ===\n\n')

global qit;

H   = gate.qft([d 1]);    % qft (generalized Hadamard) gate
add = gate.mod_add(d, d); % modular adder (generalized CNOT) gate
I   = speye(d);


% EPR preparation circuit
U = add * kron(H, I);

disp('Alice and Bob start with a shared EPR pair.')
reg = u_propagate(state('00', [d d]), U)


% two random d-its
a = floor(d*rand(1, 2));
fprintf('Alice wishes to send two d-its of information (d = %d) to Bob: a = [%d, %d].\n', d, a)


Z = diag(sqrt(d) * H(:,a(1)+1));
X = gate.mod_inc(-a(2), d);

disp('Alice encodes the d-its to her half of the EPR pair using local transformations,')
reg = u_propagate(reg, kron(Z*X, I))

disp('and sends it to Bob. He then disentangles the pair,')
reg = u_propagate(reg, U')

[p, b(1), reg] = measure(reg, 1);
[p, b(2), reg] = measure(reg, 2);
%[p, b, reg] = measure(reg, [1 2]);
b = b-1;
fprintf('and measures both qudits in the computational basis, obtaining the result  b = [%d, %d].\n', b)

if (a == b)
  disp('The d-its were transmitted succesfully.')
else
  error('Should not happen.')
end
