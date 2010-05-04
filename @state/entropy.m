function S = entropy(s, sys)
% ENTROPY  Von Neumann entropy of the state.
%  S = entropy(s)      % entropy
%  S = entropy(s ,sys) % entropy of entanglement
%
%  Returns the entropy S of the state s.
%
%  If a vector of subsystem indices sys is given, returns the
%  entropy of entanglement of the state s wrt. the partitioning
%  defined by sys.
%
%  Entropy of entanglement is only defined for pure states.
%
%  S(\rho) = -trace(\rho * \log_2(\rho))

% Ville Bergholm 2009-2010


if (nargin == 2)
  s = ptrace(to_ket(s), sys); % partial trace over one partition
end

if (size(s.data, 2) == 1)
  % pure state
  S = 0;
else
  p = eig(s.data);
  p(find(p == 0)) = 1; % avoid trouble with the logarithm
  S = -p'*log2(p);
end
