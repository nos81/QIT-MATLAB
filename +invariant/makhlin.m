function g = makhlin(c)
% MAKHLIN  Makhlin local invariants of a two-qubit gate.
%  g = makhlin(U)  % gate given as a matrix
%  g = makhlin(c)  % gate given in terms of canonical invariants
%
%  Returns a vector of three real Makhlin invariants corresponding
%  to the U(4) gate U.
%
%  Alternatively, given a row vector of canonical invariants
%  normalized to [0, 1], returns the corresponding Makhlin invariants.

%! Yu. Makhlin, "Nonlocal Properties of Two-Qubit Gates and Mixed States, and the Optimization of Quantum Computations", QIP 1, 243 (2002).
%! J. Zhang et al., "Geometric theory of nonlocal two-qubit operations", PRA 67, 042313 (2003).
% Ville Bergholm 2004-2010


if (size(c, 2) == 3)
  % matrix consisting of row vectors of canonical invariants
  c = pi*c;
  g(:,1) = (cos(c(:,1)).*cos(c(:,2)).*cos(c(:,3))).^2 -(sin(c(:,1)).*sin(c(:,2)).*sin(c(:,3))).^2;
  g(:,2) = 0.25*sin(2*c(:,1)).*sin(2*c(:,2)).*sin(2*c(:,3));
  g(:,3) = 4*g(:,1) - cos(2*c(:,1)).*cos(2*c(:,2)).*cos(2*c(:,3));
else
  % gate matrix    
  U = c;
  Q_Bell = [1 0 0 i; 0 i 1 0; 0 i -1 0; 1 0 0 -i] / sqrt(2);
  V = Q_Bell' * U * Q_Bell;
  M = V.'*V;

  t1 = trace(M)^2;
  t2 = t1 / (16*det(U));
  g(1) = real(t2);
  g(2) = imag(t2);
  g(3) = real((t1 - trace(M*M)) / (4*det(U)));
end
