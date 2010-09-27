function F = fidelity(r, s)
% FIDELITY  Fidelity of two states.
%  F = fidelity(rho, sigma)
%
%  Fidelity of two state operators \rho and \sigma is defined as
%  $F(\rho, \sigma) = \trace \sqrt{\sqrt{\rho} \sigma \sqrt{\rho}}$.
%  For state vectors this is equivalent to the overlap, F = |<a|b>|.
%
%  Fidelity is symmetric in its arguments and bounded in the interval [0,1].

% Ville Bergholm 2009-2010
%! M.A. Nielsen, I.L. Chuang, "Quantum Computation and Quantum Information" (2000), chapter 9.2.2


if (size(r.data, 2) == 1)
  if (size(s.data, 2) == 1)
    F = abs(r.data' * s.data);
  else
    F = sqrt(real(r.data' * s.data * r.data));
  end
else
  if (size(s.data, 2) == 1)
    F = sqrt(real(s.data' * r.data * s.data));
  else
    temp = sqrtm(r.data);
    F = real(trace(sqrtm(temp * s.data * temp)));
  end
end
