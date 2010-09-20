function s = reorder(s, p)
% REORDER  Change the order of subsystems in a state.
%  x = reorder(s, perm);
%  x = reorder(s, [2 4 3 1]); % reorder the subsystems of s to the order given
%  x = reorder(s, [2 5]);     % swap subsystems 2 and 5
%
%  Reorders the subsystems of state s according to permutation vector perm.
%  If only two subsystems are listed, swaps them.

% Ville Bergholm 2009


if (nargin < 2)
  error('Need the permutation.')
end

% this is just a wrapper for utils/reorder
[s.data, s.dim] = reorder(s.data, s.dim, p);
