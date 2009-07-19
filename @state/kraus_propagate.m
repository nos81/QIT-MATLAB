function s = kraus_propagate(s, E)
% STATE/KRAUS_PROPAGATE  Apply a quantum operation on a state.
%  s1 = kraus_propagate(s0, E)
%
%  Applies the quantum operation E on the state s0.
%  E == {E_1, E_2, ...} is a set of Kraus operators.

% Ville Bergholm 2009


if (nargin < 2)
  error('Needs the quantum operation.');
end


if (size(s.data, 1) ~= size(E{1}, 2))
  error('Dimension of the operation does not match the dimension of the state.');
  % TODO allow the user to apply E only to some subsystems of s0
end


n = length(E);
% TODO: If n > prod(s.dim)^2, there is a simpler equivalent
% operation. Should the user be notified?

temp = 0;
for k=1:n
  temp = temp + E{k}'*E{k};
end
% TODO inform user if E is unphysical?


if (size(s.data, 2) == 1)
  % state vector
  if (n == 1)
    s.data = E{1}*s.data; % remains a pure state
    return;
  end
  s.data = s.data*s.data'; % into a state operator
end

% state operator
temp = 0;
for k=1:n
  temp = temp + E{k}*s.data*E{k}';
end
s.data = temp;
