function p = phase_estimation_precision(t, U, u)
% PHASE_ESTIMATION_PRECISION  Quantum phase estimation demo.
%
%  p = phase_estimation_precision(t, U [, u])
%
%  Estimate an eigenvalue of unitary operator U using t bits, starting from the state u.
%  Plots and returns the probability distribution of the resulting t-bit approximations.
%  If u is not given, use a random eigenvector of U.

%! R. Cleve et al., "Quantum Algorithms Revisited", Proc. R. Soc. London A454, 339 (1998).
%! M.A. Nielsen, I.L. Chuang, "Quantum Computation and Quantum Information" (2000), chapter 5.2.
% Ville Bergholm 2009-2010


fprintf('\n\n=== Phase estimation ===\n\n')

% find eigenstates of the operator
if (nargin < 2)
  error('Both t and U must be given.')
end

N = size(U, 1);
[v,d] = eig(U);
if (nargin < 3)
  u = state(v(:,1), N); % exact eigenstate
end

fprintf('Use %d qubits to estimate the phases of the eigenvalues of a U(%d) operator.\n', t, N)

p = real(prob(phase_estimation(t, U, u))); % fix rounding errors

T = 2^t;

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
