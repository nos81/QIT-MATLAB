function p = purify(s)
% STATE/PURIFY  Purifies the state, converts the representation into a ket.
%  p = purify(s)
%
%  TODO: At present this function does NOT create a purification of an arbitrary mixed state s.
%  It simply converts a state operator into the corresponding state vector if the state is pure.

% Ville Bergholm 2009


global qit;

p = s;

if (size(s.data, 2) == 1)
  return; % nothing to do
end

if (abs(purity(s) - 1) > qit.tol)
  error('The state is not pure, and cannot be represented using a ket vector.')
  return;
end

[v, d] = eig(p.data);
v = v(:,end); % corresponds to the highest eigenvalue, e.g. 1

% phase convention: first nonzero element in state vector is real, positive
for k=1:length(v)
  if (abs(v(k)) > qit.tol)
    phase = v(k)/abs(v(k));
    v = v / phase;
    break;
  end
end

p.data = v;
