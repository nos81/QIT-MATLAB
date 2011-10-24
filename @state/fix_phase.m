function s = fix_phase(s)
% FIX_PHASE  Apply a global phase convention to a ket state.
%
%  q = fix_phase(s)
%
%  If the state s is represented using a ket, returns a copy of s for
%  which the global phase has been set such that the first nonzero
%  element in the state vector is real and positive.

% Ville Bergholm 2009


global qit;

if ~is_ket(s)
  % state operator
  return;
end

% apply the phase convention: first nonzero element in state vector is real, positive
v = s.data;
for k=1:length(v)
  if (abs(v(k)) > qit.tol)
    phase = v(k)/abs(v(k));
    s.data = v / phase;
    return;
  end
end
