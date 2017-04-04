function bb84(n)
% BB84  Bennett-Brassard 1984 quantum key distribution protocol demo.
%
%  bb84(n)
%
%  Simulate the protocol with n qubits transferred.

% Ville Bergholm 2009-2016


fprintf('\n\n=== BB84 protocol ===\n\n')

global qit;

if nargin < 1
  n = 50; % string len
end

fprintf('Using %d transmitted qubits, intercept-resend attack.\n\n', n)

H  = qit.H;  % Hadamard gate
sx = qit.sx; % bit flip

% Alice generates two random bit vectors
bits_A  = rand(1, n) > 0.5;
basis_A = rand(1, n) > 0.5;

% Bob generates one random bit vector
basis_B = rand(1, n) > 0.5;

fprintf(['Alice transmits a sequence of qubits to Bob using a quantum channel.\n',...
         'For every qubit, she randomly chooses a basis (computational or diagonal)\n',...
         'and randomly prepares the qubit in either the |0> or the |1> state in that basis.\n'])
fprintf('\nAlice''s bits: \t')
fprintf('%c', bits_A+'0')
fprintf('\nAlice''s basis: \t')
fprintf('%c', basis_A+'0')

temp = state('0');
for k=1:n
  % Alice has a source of qubits in the zero state.
  q = temp;

  % Should Alice flip the qubit?
  if bits_A(k), q = q.prop(sx); end

  % Should Alice apply a Hadamard?
  if basis_A(k), q = q.prop(H); end

  % Alice now sends the qubit to Bob...
  % ===============================================
  % ...but Eve intercepts it, and conducts an intercept/resend attack

  % Eve might have different strategies here... TODO
  % Eve's strategy (simplest choice, random)
  basis_E(k) = rand(1) > 0.5;
  
  % Eve's choice of basis: Hadamard?
  if basis_E(k), q = q.prop(H); end

  % Eve measures in the basis she has chosen
  [p, res, q] = measure(q);
  bits_E(1, k) = res-1;

  % Eve tries to reverse the changes she made...
  if basis_E(k), q = q.prop(H); end

  % ...and sends the result to Bob.
  % ===============================================

  % Bob's choice of basis
  if basis_B(k), q = q.prop(H); end

  % Bob measures in the basis he has chosen, and discards the qubit.
  [p, res] = measure(q);
  bits_B(1, k) = res-1;
end


fprintf(['\n\n\nHowever, there''s an evesdropper, Eve, on the line. She intercepts the qubits,\n',...
         'randomly measures them in either basis (thus destroying the originals!), and then sends\n',...
         'a new batch of qubits corresponding to her measurements and basis choices to Bob.\n',...
         'Since Eve on the average can choose the right basis only 1/2 of the time,\n',...
         'about 1/4 of her bits differ from Alice''s.\n'])

fprintf('\nEve''s basis: \t')
fprintf('%c', basis_E+'0')
fprintf('\nEve''s bits: \t')
fprintf('%c', bits_E+'0')
fprintf('\n\nMismatch frequency between Alice and Eve: %.3f\n', sum(xor(bits_A, bits_E))/n)

fprintf('\n\nWhen Bob receives the qubits, he randomly measures them in either basis.\n');
fprintf('Due to Eve''s eavesdropping, Bob''s bits differ from Alice''s 3/8 of the time.\n');
fprintf('\nBob''s basis: \t')
fprintf('%c', basis_B+'0')
fprintf('\nBob''s bits: \t')
fprintf('%c', bits_B+'0')
fprintf('\n\nMismatch frequency between Alice and Bob: %.3f\n', sum(xor(bits_A, bits_B))/n)

%sum(xor(bits_A, bits_E))/n
%sum(xor(bits_A, bits_B))/n

fprintf(['\n\nNow Bob announces on a public classical channel that he has received all the qubits.\n',...
         'Alice and Bob then reveal the bases they used. Whenever the bases happen to match,\n',...
         '(about 1/2 of the time on the average), they both add their corresponding bit to\n'...
         'their personal key. The two keys should be identical unless there''s been an evesdropper.\n'...
         'However, because of Eve each key bit has a 1/4 probability of being wrong.\n']);

match = not(xor(basis_A, basis_B));

fprintf('\n%d basis matches.', sum(match))
key_A = bits_A(find(match));
key_B = bits_B(find(match));

fprintf('\nAlice''s key: \t')
fprintf('%c', key_A+'0')
fprintf('\nBob''s key: \t')
fprintf('%c', key_B+'0')

m = length(key_A);
fprintf('\n\nMismatch frequency between Alice''s and Bob''s keys: %.3f\n\n', sum(xor(key_A, key_B))/m)

fprintf(['Alice and Bob then sacrifice k bits of their shared key to compare them publicly.\n',...
         'If a nonmatching bit is found, the reason is either an eavesdropper or a noisy channel.\n',...
         'With a noiseless channel Alice and Bob will detect Eve''s presence with the probability 1-(3/4)^k.\n'])
