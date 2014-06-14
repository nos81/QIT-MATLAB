function S = entropy(s, sys, alpha)
% ENTROPY  Von Neumann or Renyi entropy of the state.
%  S = entropy(s [, sys=[] [, alpha=1]])
%
%  Returns the Renyi entropy of order alpha,
%  S_\alpha(\rho) = \frac{1}{1-\alpha} \log_2 \trace(\rho^\alpha).
%
%  When \alpha == 1 (default), this coincides with the von Neumann entropy
%  S(\rho) = -\trace(\rho \log_2(\rho)).
%
%  If sys == [], returns the entropy of the state s.
%  If sys is a vector of subsystem indices, returns the
%  entropy of entanglement of the state s wrt. the partitioning
%  defined by sys. Entropy of entanglement is only defined for pure states.

% Ville Bergholm 2009-2014

if nargin < 3
    alpha = 1;
end
if nargin < 2
    sys = [];
end

if ~isempty(sys)
  s = ptrace(to_ket(s), sys); % partial trace over one partition
end

if is_ket(s)
  % pure state
  S = 0;
else
  p = eig(s.data);
  if alpha ~= 1
      % Renyi entropy
      S = log2(sum(p .^ alpha)) / (1-alpha);
  else
      % Von Neumann entropy
      p(find(p == 0)) = 1; % avoid trouble with the logarithm
      S = -p'*log2(p);
  end
end
