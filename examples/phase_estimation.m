function p = phase_estimation(t, U, u)
% PHASE_ESTIMATION  Quantum phase estimation algorithm.
%
%  p = phase_estimation(t, U [, u])
%
%  Estimate an eigenvalue of unitary operator U using t bits, starting from the state u.
%  Returns the probability distribution of the resulting t-bit approximations.
%  If u is not given, use a random eigenvector of U.

%! R. Cleve et al., "Quantum Algorithms Revisited", Proc. R. Soc. London A454, 339 (1998).
%! M.A. Nielsen, I.L. Chuang, "Quantum Computation and Quantum Information" (2000), chapter 5.2.
% Ville Bergholm 2009


% TODO to get a result accurate to n bits with probability >= (1-epsilon),
% choose  t >= n + ceil(log2(2+1/(2*epsilon)))

fprintf('\n\n=== Phase estimation ===\n\n')

% find eigenstates of the operator
if (nargin < 2)
  error('Both t and U must be given.')
end

[v,d] = eig(U);
if (nargin < 3)
  u = v(:,1); % exact eigenstate
end

T = 2^t;
N = size(U, 1);

fprintf('Use %d qubits to estimate the phases of the eigenvalues of a U(%d) operator.\n', t, N)

% index register
reg_t = state(ones(T, 1)/sqrt(T), 2*ones(1,t)); % uniform superposition

% state register
reg_u = state(u, N);

% full register
reg = tensor(reg_t, reg_u);

% controlled unitaries
for k = 1:t
    ctrl = -ones(1, t);
    ctrl(k) = 1;
    temp = gate.controlled(U^(2^(t-k)), ctrl);
    
    reg.data = temp * reg.data;
end

% inverse QFT
reg.data = kron(gate.qft(t)', eye(N)) * reg.data;

% measurements
reg_t = ptrace(reg, t+1); % trace over state register
p = prob(reg_t); % measurement probabilities

% TODO lighter formula
% [p] = measure(reg, 1:t)
%xxx = reg.data;
%for k = 1:T
%  temp = N*(k-1);  
%  p(k) = norm(xxx(temp+1:temp+N))^2;
%end

% plot probability distribution
figure;
hold on;
bar((0:T-1)/T, p);
xlabel('phase/2\pi');
ylabel('probability');
title('Phase estimation')
axis([-1/(T*2) 1-1/(T*2) 0 1])

% compare to correct answer
target = angle(diag(d))/(2*pi) + 1;
target = target - floor(target)
plot(target, 0.5*max(p)*ones(size(target)), 'mo');

legend('Measurement probability distribution', 'Target phases');
