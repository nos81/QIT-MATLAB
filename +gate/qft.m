function U = qft(dim)
% QFT  Quantum Fourier transform gate.
%  U = qft(dim)
%
%  Returns the quantum Fourier transform gate for the specified system.
%  dim is a vector of subsystem dimensions.

% Ville Bergholm 2004-2010


global qit;

n = length(dim);
cache = false;

if (isequal(dim, qubits(n)))
  cache = true;

  % check cache first
  if (length(qit.qft) >= n && length(qit.qft{n}) > 0)
    U = qit.qft{n};
    return;
  end
end

N = prod(dim);

U = zeros(N);
for j=1:N
  for k=1:N
    U(j,k) = exp(i*2*pi*(j-1)*(k-1)/N) / sqrt(N);
  end
end

U = lmap(U, {dim, dim});

if (cache)
  % store it in the cache
  qit.qft{n} = U;
end
