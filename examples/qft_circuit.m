function U = qft_circuit(dim)
% QFT_CIRCUIT  Quantum Fourier transform circuit demo.
%
%  U = qft_circuit(dim)
%
%  Simulate the quadratic QFT circuit construction.
%  dim is the dimension vector of the subsystems.
%
%  NOTE: If dim is not palindromic the resulting circuit also
%  reverses the order of the dimensions

% U |x1,x2,...,xn> = 1/sqrt(d) \sum_{ki} |kn,...,k2,k1> exp(i 2 \pi (k1*0.x1x2x3 +k2*0.x2x3 +k3*0.x3))
% = 1/sqrt(d) \sum_{ki} |kn,...,k2,k1> exp(i 2 \pi 0.x1x2x3*(k1 +d1*k2 +d1*d2*k3))

% Ville Bergholm 2010


fprintf('\n\n=== Quantum Fourier transform using a quadratic circuit ===\n\n')

if (nargin < 1)
  dim = [2 3 3 2]
end

n = length(dim)
U = gate.id(dim);

for k=1:n
  H = gate.qft(dim(k));
  U = gate.single(H, k, dim) * U;
  for j=k+1:n
    temp = Rgate(dim(k:j));
    U = gate.two(temp, [k, j], dim) * U;
  end
end

for k=1:n/2
  temp = gate.swap(dim(k), dim(n+1-k));
  U = gate.two(temp, [k, n+1-k], dim) * U;
end

U.data = full(U.data); % it's a QFT anyway

U
%temp = U - gate.qft(dim);
%norm(temp.data)
end


function R = Rgate(d)
% R = \sum_{xy} exp(i*2*pi * x*y/prod(dim)) |xy><xy|

  temp = kron(0:d(1)-1, 0:d(end)-1)/prod(d);
  R = gate.phase(2*pi*temp, [d(1), d(end)]);
end
