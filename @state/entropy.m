function S = entropy(s)
% ENTROPY  Von Neumann entropy of the state.
%  S = entropy(s)
%
%  S(\rho) = -trace(\rho * \log_2(\rho))

% Ville Bergholm 2009


if (size(s.data, 2) == 1)
  % pure state
  S = 0;
else
  p = eig(s.data);
  p(find(p == 0)) = 1; % avoid trouble with the logarithm
  S = -p'*log2(p);
end
