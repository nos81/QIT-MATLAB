function s = to_ket(s)
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
  return; % nothing to do
end

if (abs(purity(s) - 1) > qit.tol)
  error('The state is not pure, and thus cannot be represented using a ket vector.')
end

[v, d] = eig(s.data);
v = v(:,end); % corresponds to the highest eigenvalue, e.g. 1

% phase convention: first nonzero element in state vector is real, positive
for k=1:length(v)
  if (abs(v(k)) > qit.tol)
    phase = v(k)/abs(v(k));
    v = v / phase;
    break;
  end
end

s.data = v;
