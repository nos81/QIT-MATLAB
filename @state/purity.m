function p = purity(s)
% STATE/PURITY  Purity of the state.
%  p = purity(s)
%
%  Returns the purity of the state, p = trace(\rho^2).
%  Equivalent to linear entropy, S_l = 1-p.
%  The state is first normalized. TODO reasonable?

% Ville Bergholm 2008-2009


if (size(s.data, 2) == 1)
  p = 1;
else
  temp = s.data / trace(s.data); % normalize it
  p = trace(temp * temp);
end
