function s = to_ket(s, set_phase)
% STATE/TO_KET  Convert state representation into a ket (if possible).
%
%  q = to_ket(s)
%
%  If the state s is pure, returns q, a copy of s for which the
%  internal representation of the state (q.data) is guaranteed to
%  be a ket vector.

% Ville Bergholm 2009


global qit;

if (size(s.data, 2) == 1)
  % already a ket
  if (nargin < 2)
    return; % nothing to do
  else
    v = s.data; % prepare to apply phase convention
  end
else
  % state op
  if (abs(purity(s) - 1) > qit.tol)
    error('The state is not pure, and thus cannot be represented using a ket vector.')
  end

  [v, d] = eig(s.data);
  v = v(:,end); % corresponds to the highest eigenvalue, i.e. 1
end

% apply the phase convention: first nonzero element in state vector is real, positive
for k=1:length(v)
  if (abs(v(k)) > qit.tol)
    phase = v(k)/abs(v(k));
    s.data = v / phase;
    return;
  end
end
