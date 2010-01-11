function s = normalize(s)
% STATE/NORMALIZE  Normalize the state.
%  q = normalize(s)
%
%  Returns the state s normalized.

% Ville Bergholm 2008

if (size(s.data, 2) == 1)
  s.data = s.data/norm(s.data);
else
  s.data = s.data/trace(s.data);
end
