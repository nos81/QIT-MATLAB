function [g1, g2, g3] = makhlin(c1, c2, c3)
% MAKHLIN  Compute the Makhlin gate invariants of a two-qubit gate.
%  [g1, g2, g3] = makhlin(U)  % gate given as a matrix
%  [g1, g2, g3] = makhlin(c1, c2, c3)  % gate given in terms of canonical invariants
%
%  Returns a vector of three real Makhlin invariants for the U(4) gate U.

%! Yu. Makhlin, "Nonlocal Properties of Two-Qubit Gates and Mixed States, and the Optimization of Quantum Computations", QIP 1, 243 (2002).
%! J. Zhang et al., "Geometric theory of nonlocal two-qubit operations", PRA 67, 042313 (2003).
% Ville Bergholm 2004-2008


if (nargin == 3)
  % canonical invariants
  s = size(c1);
  g1 = (cos(c1).*cos(c2).*cos(c3)).^2 -(sin(c1).*sin(c2).*sin(c3)).^2;
  g2 = 0.25*sin(2*c1).*sin(2*c2).*sin(2*c3);
  g3 = 4*g1 - cos(2*c1).*cos(2*c2).*cos(2*c3);
else
  % gate matrix    
  U = c1;
  Q_Bell = [1 0 0 i; 0 i 1 0; 0 i -1 0; 1 0 0 -i] / sqrt(2);
  V = Q_Bell' * U * Q_Bell;
  M = V.'*V;

  t1 = trace(M)^2;
  t2 = t1 / (16*det(U));
  g1 = real(t2);
  g2 = imag(t2);
  g3 = (t1 - trace(M*M)) / (4*det(U));
end
