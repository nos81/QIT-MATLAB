function [s] = reorder(s, perm)
% REORDER  Change the relative order of subsystems in a state.
%  q = reorder(s, perm);
%    reorder(s, [3 2 1]); % reverse the order of subsystems
%    reorder(s, [2 5]);   % swap subsystems 2 and 5
%
%  Reorders the subsystems of the state s according to permutation vector perm.
%
%  The permutation vector may consist of either exactly two subsystem numbers
%  (to be swapped), or a full permutation of subsystem numbers.

% Ville Bergholm 2010


% this is just an adapter for lmap::reorder

if (size(s.data, 2) == 1)
  s = reorder@lmap(s, {perm});
else
  s = reorder@lmap(s, {perm, perm});
end
