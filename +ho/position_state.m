function s = position_state(x, n)
% POSITION_STATE  Position eigenstates of a harmonic oscillator.
%  s = position_state(x, n)
%
%  Returns the n-dimensional approximation of the eigenstate |x>
%  of the dimensionless position operator X in the number basis.
%  X = \sqrt{\frac{2 m \omega}{\hbar}} x = a+a'.
%
%  H = \frac{p^2}{2m} +\frac{1}{2} m \omega^2 x^2 = \frac{1}{4} \hbar \omega (P^2 +X^2)

% Ville Bergholm 2010


ket = zeros(n, 1);

ket(1) = 1; % arbitrary
ket(2) = x*ket(1);
for k=3:n
  n = k-1; % occupation number for ket(k) (MATLAB indexing!)
  ket(k) = x/sqrt(n) * ket(k-1) -sqrt((n-1)/n) * ket(k-2);
end

ket = ket/norm(ket);
s = state(ket, n);
