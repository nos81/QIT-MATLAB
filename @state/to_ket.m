function s = to_ket(s)
% TO_KET  Convert state representation into a ket (if possible).
%
%  q = to_ket(s)
%
%  If the state s is pure, returns a copy of s for which the
%  internal representation of the state (q.data) is guaranteed to
%  be a ket vector.

% Ville Bergholm 2009


global qit;

if (size(s.data, 2) == 1)
  return; % already a ket, nothing to do
else
  % state op
  if (abs(purity(s) - 1) > qit.tol)
    error('The state is not pure, and thus cannot be represented using a ket vector.')
  end

  [v, d] = eig(s.data);
  s.data = v(:,end); % corresponds to the highest eigenvalue, i.e. 1
  s = fix_phase(s); % clean up global phase
end
