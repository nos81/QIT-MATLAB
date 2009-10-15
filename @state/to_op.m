function s = to_op(s)
% STATE/TO_OP  Use the state operator representation for a quantum state.
%
%  q = to_op(s)
%
%  Returns q, a copy of s for which the internal representation 
%  of the state (q.data) is guaranteed to be a state operator.

% Ville Bergholm 2009


if (size(s.data, 2) ~= 1)
  return; % nothing to do
end

s.data = s.data*s.data';
