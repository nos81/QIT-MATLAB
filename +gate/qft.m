function U = qft(n)
% QFT  Quantum Fourier transform gate.
%  U = qft(n)
%
%  Returns the quantum Fourier transform matrix for the specified system.
%  If n is a scalar, the system consists of n qubits.
%  Otherwise n is assumed to be a vector of subsystem dimensions.
%  A single d-dimensional system can be specified using n = [d 1].

% Ville Bergholm 2004-2009

global qit;

if (isscalar(n))
  N = 2^n; % n qubits

  % check cache first
  if (length(qit.qft) >= n && length(qit.qft{n}) > 0)
    U = qit.qft{n};
    return;
  end

else
  N = prod(n); % vector of dimensions
end


U = zeros(N);
for j=1:N
  for k=1:N
    U(j,k) = exp(i*2*pi*(j-1)*(k-1)/N) / sqrt(N);
  end
end

if (isscalar(n))
  % store it in the cache
  qit.qft{n} = U;
end
