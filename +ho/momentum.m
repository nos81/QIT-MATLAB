function P = momentum(n)
% MOMENTUM  Momentum operator.
%  P = momentum(n)
%
%  Returns the n-dimensional approximation of the
%  dimensionless momentum operator P in the number basis.
%
%  Q = \sqrt{\frac{2 m \omega}{\hbar}} q =     a+a'
%  P = \sqrt{\frac{2}{m \hbar \omega}} p = -i*(a-a')
%
%  [q, p] = i \hbar,  [Q, P] = 2i
%
%  H = \frac{p^2}{2m} +\frac{1}{2} m \omega^2 q^2
%    = 1/4 \hbar \omega (P^2 +Q^2)
%    = \hbar \omega (a'*a +\frac{1}{2})

% Ville Bergholm 2010


a = boson_ladder(n);
P = -1i*(a-a');
