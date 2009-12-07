function bb84(n)
% BB84  Bennett-Brassard 1984 quantum key distribution protocol demo.
%
%  bb84(n)
%
%  Simulate the protocol with n qubits transferred.

% Ville Bergholm 2009


fprintf('\n\n=== BB84 protocol ===\n\n')

global qit;

if (nargin < 1)
  n = 50; % string len
end

fprintf('Using %d transmitted qubits.\n\n', n)

H  = qit.H;  % Hadamard gate
sx = qit.sx; % bit flip

% Alice generates two random bit vectors
sent    = rand(1, n) > 0.5;
basis_A = rand(1, n) > 0.5;

% Bob generates one random bit vector
basis_B = rand(1, n) > 0.5;

disp('Alice transmits a sequence of qubits to Bob using a quantum channel.')
disp('For every qubit, she randomly chooses a basis (computational or diagonal)')
disp('and randomly prepares the qubit in either the |0> or the |1> state in that basis.')
fprintf('\n')
disp('When Bob receives the qubits, he randomly measures them in either basis');

temp = state('0');
for k=1:n
  % Alice has a source of qubits in the zero state.
  q = temp;

  % Should Alice flip the qubit?
  if (sent(k)), q = u_propagate(q, sx); end

  % Should Alice apply a Hadamard?
  if (basis_A(k)), q = u_propagate(q, H); end

  % Alice now sends the qubit to Bob...
  % ===============================================
  % ...but Eve intercepts it, and conducts an intercept/resend attack

  % Eve might have different strategies here... TODO
  % Eve's strategy (simplest choice, random)
  basis_E(k) = rand(1) > 0.5;
  
  % Eve's choice of basis: Hadamard?
  if (basis_E(k)), q = u_propagate(q, H); end

  % Eve measures in the basis she has chosen
  [p, res, q] = measure(q);
  eavesdrop(1, k) = res-1;

  % Eve tries to reverse the changes she made...
  if (basis_E(k)), q = u_propagate(q, H); end

  % ...and sends the result to Bob.
  % ===============================================

  % Bob's choice of basis
  if (basis_B(k)), q = u_propagate(q, H); end

  % Bob measures in the basis he has chosen, and discards the qubit.
  [p, res] = measure(q);
  received(1, k) = res-1;
end

%sum(xor(sent, eavesdrop))/n
%sum(xor(sent, received))/n

disp('Now Bob announces on a public classical channel that he has received all the qubits.')
disp('Alice then reveals the bases she used, and Bob compares them to his.')
disp('Whenever the bases match, so should the prepared/measured values unless there''s an eavesdropper.')

match = not(xor(basis_A, basis_B));

key_A = sent(find(match))
key_B = received(find(match))
m = length(key_A);
fprintf('\nMismatch frequency between Alice and Bob: %f\n\n', sum(xor(key_A, key_B))/m)

disp('Alice and Bob then sacrifice k bits of their shared key to compare them.')
disp('If an nonmatching bit is found, the reason is either an eavesdropper or a noisy channel.')
disp('Since the probability for each eavesdropped bit to be wrong is 1/4, they will detect.')
disp('Eve''s presence with the probability 1-(3/4)^k.')
