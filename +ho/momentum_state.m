function s = momentum_state(p, n)
% MOMENTUM_STATE  Momentum eigenstates of a harmonic oscillator.
%  s = momentum_state(p, n)
%
%  Returns the n-dimensional approximation of the eigenstate |p>
%  of the dimensionless momentum operator P in the number basis.
%  P = \sqrt{\frac{2}{\hbar m \omega}} p = i*(a-a').
%
%  H = \frac{p^2}{2m} +\frac{1}{2} m \omega^2 x^2 = \frac{1}{4} \hbar \omega (P^2 +X^2)

% Ville Bergholm 2010


ket = zeros(n, 1);

ket(1) = 1; % arbitrary
ket(2) = i*x*ket(1);
for k=3:n
  n = k-1; % occupation number for ket(k) (MATLAB indexing!)
  ket(k) = i*x/sqrt(n) * ket(k-1) +sqrt((n-1)/n) * ket(k-2);
end

ket = ket/norm(ket);
s = state(ket, n);
