function q = kraus_propagate(s, E)
% KRAUS_PROPAGATE  Apply a quantum operation to the state.
%  s1 = kraus_propagate(s0, E)
%
%  Applies the quantum operation E on the state s0.
%  E == {E_1, E_2, ...} is a set of Kraus operators.

% Ville Bergholm 2009


global qit;

if (nargin < 2)
  error('Needs the quantum operation.');
end


% TODO allow the user to apply E only to some subsystems of s0

n = length(E);
% TODO: If n > prod(dims(s))^2, there is a simpler equivalent
% operation. Should the user be notified?

%temp = 0;
%for k=1:n
%  temp = temp + E{k}'*E{k};
%end
%if (norm(temp - eye(size(temp))) > qit.tol)
%  warning('Unphysical quantum operation.')
%end

if is_ket(s)
  % state vector
  if (n == 1)
    q = u_propagate(s, E{1}); % remains a pure state
    return;
  end
  s = to_op(s); % into a state operator
end

% state operator
q = u_propagate(s, E{1});
for k=2:n
  q = q + u_propagate(s, E{k});
end
