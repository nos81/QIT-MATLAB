function s = normalize(s)
% NORMALIZE  Normalize the state.
%  q = normalize(s)
%
%  Returns the state s normalized.

% Ville Bergholm 2008

if is_ket(s)
  s.data = s.data/norm(s.data);
else
  s.data = s.data/trace(s.data);
end
