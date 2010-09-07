function p = purity(s)
% PURITY  Purity of the state.
%  p = purity(s)
%
%  Returns the purity of the state, p = trace(\rho^2).
%  Equivalent to linear entropy, S_l = 1-p.

% Ville Bergholm 2008-2010


if (size(s.data, 2) == 1)
  p = 1;
else
  p = trace(s.data^2);
end
