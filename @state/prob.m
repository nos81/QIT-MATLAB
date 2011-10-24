function p = prob(s)
% PROB  Measurement probabilities of the state in the computational basis.
%  p = prob(s)
%
%  Vector p holds the probabilities of finding a system
%  in the state s in all the different states of the computational basis.

% Ville Bergholm 2009


if is_ket(s)
  % state vector
  p = abs(s.data).^2;
else
  % state operator
  p = diag(s.data);
end
