function s = momentum_state(p, n)
% MOMENTUM_STATE  Momentum eigenstates of a harmonic oscillator.
%  s = momentum_state(p, n)
%
%  Returns the n-dimensional approximation of the eigenstate |p>
%  of the dimensionless momentum operator P in the number basis.
%
%  See position.m, momentum.m.

% Ville Bergholm 2010


ket = zeros(n, 1);

temp = 1i*p;

ket(1) = 1; % arbitrary
ket(2) = temp * ket(1);
for k=3:n
  nn = k-1; % occupation number for ket(k) (MATLAB indexing!)
  ket(k) = temp/sqrt(nn) * ket(k-1) +sqrt((nn-1)/nn) * ket(k-2);
end

ket = ket/norm(ket);
s = state(ket, n);
