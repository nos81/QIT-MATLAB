function p = prob(s)
% STATE/PROB  Measurement probabilities of the state in the computational basis.
%  p = prob(s)
%
%  Vector p holds the probabilities of finding a system
%  in the state s in all the different states of the computational basis.

% Ville Bergholm 2009


if (size(s.data, 2) == 1)
  % state vector
  p = abs(s.data).^2;
else
  % state operator
  p = diag(s.data);
end
