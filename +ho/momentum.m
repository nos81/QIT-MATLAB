function P = momentum(n)
% MOMENTUM  Momentum operator.
%  P = momentum(n)
%
%  Returns the n-dimensional approximation of the
%  dimensionless momentum operator P in the number basis.
%
%  Q = \sqrt{\frac{m \omega}{\hbar}}   q =    (a+a')/sqrt(2)
%  P = \sqrt{\frac{1}{m \hbar \omega}} p = -i*(a-a')/sqrt(2)
%
%  [q, p] = i \hbar,  [Q, P] = i
%
%  H = \frac{p^2}{2m} +\frac{1}{2} m \omega^2 q^2
%    = 0.5 \hbar \omega (P^2 +Q^2)
%    = \hbar \omega (a'*a +\frac{1}{2})

% Ville Bergholm 2010


a = boson_ladder(n);
P = -i*(a-a')/sqrt(2);
