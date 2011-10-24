function s = to_op(s)
% TO_OP  Convert state representation into a state operator.
%
%  q = to_op(s)
%
%  Returns q, a copy of s for which the internal representation 
%  of the state (q.data) is guaranteed to be a state operator.

% Ville Bergholm 2009-2011


if ~is_ket(s)
  return; % nothing to do
end

s.data = s.data * s.data';
s.dim{2} = s.dim{1};
