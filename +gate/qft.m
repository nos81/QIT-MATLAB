function U = qft(n)
% QFT  Quantum Fourier transform gate.
%  U = qft(n)
%
%  Returns the quantum Fourier transform matrix for n qubits.

% Ville Bergholm 2004-2009

N = 2^n;
U = zeros(N);

for j=1:N
  for k=1:N
    U(j,k) = exp(i*2*pi*(j-1)*(k-1)/N) / sqrt(N);
  end
end
