function reg = phase_estimation(t, U, reg_s, implicit)
% PHASE_ESTIMATION  Quantum phase estimation algorithm.
%  reg = phase_estimation(t, U, s)
%
%  Estimate an eigenvalue of unitary operator U using t qubits,
%  starting from the state s.
%
%  Returns the state of the index register after the phase estimation
%  circuit, but before final measurement.
%
%  To get a result accurate to n bits with probability >= (1-epsilon),
%  choose  t >= n + ceil(log2(2+1/(2*epsilon))).

%! R. Cleve et al., "Quantum Algorithms Revisited", Proc. R. Soc. London A454, 339 (1998).
%! M.A. Nielsen, I.L. Chuang, "Quantum Computation and Quantum Information" (2000), chapter 5.2.
% Ville Bergholm 2009-2010


% find eigenstates of the operator
if (nargin < 3)
  error('Need t, U and s.')
end

T = 2^t;
S = size(U, 1);

% index register in uniform superposition
%reg_t = u_propagate(state(0, qubits(t)), gate.walsh(t)); % use Hadamards
reg_t = state(ones(T, 1)/sqrt(T), qubits(t)); % skip the Hadamards

% state register (ignore the dimensions)
reg_s = state(reg_s, S);

% full register
reg = tensor(reg_t, reg_s);

% controlled unitaries
for k = 1:t
    ctrl = -ones(1, t);
    ctrl(k) = 1;
    temp = gate.controlled(U^(2^(t-k)), ctrl);
    reg = u_propagate(reg, temp);
end

% from this point forward the state register is not used anymore

if (nargin == 4 && implicit)
  % save memory and CPU: make an implicit measurement of the state reg, discard the results
  [dummy, res, reg] = measure(reg, t+1, true);
  %fprintf('Implicit measurement of state register: %d\n', res);
else
  % more expensive computationally: trace over the state register
  reg = ptrace(reg, t+1);
end

% do an inverse quantum Fourier transform on the index reg
QFT = gate.qft(qubits(t));
reg = u_propagate(reg, QFT');
