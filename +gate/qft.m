function U = qft(n)
% QFT  Quantum Fourier transform gate.
%  U = qft(n)
%
%  Returns the quantum Fourier transform matrix for the specified system.
%  If n is a scalar, the system consists of n qubits.
%  Otherwise n is assumed to be a vector of subsystem dimensions.
%  A single d-dimensional system can be specified using n = [d 1].

% Ville Bergholm 2004-2009


if (isscalar(n))
  N = 2^n; % n qubits
else
  N = prod(n); % vector of dimensions
end

U = zeros(N);
for j=1:N
  for k=1:N
    U(j,k) = exp(i*2*pi*(j-1)*(k-1)/N) / sqrt(N);
  end
end
