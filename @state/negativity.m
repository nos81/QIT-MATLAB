function N = negativity(s, sys)
% NEGATIVITY  Negativity of the state.
%  N = negativity(s, sys)
%
%  Returns the negativity of the state s wrt. the partitioning
%  given by the listing of subsystems in the vector sys.

%! A. Peres, "Separability Criterion for Density Matrices", PRL 77, 1413 (1996).
%! M. Horodecki et al., "Separability of Mixed States: Necessary and Sufficient Conditions", Physics Letters A 223, 1-8 (1996).
% Ville Bergholm 2008-2014


if (nargin < 2)
  error('Need the partitioning.')
end

s = ptranspose(s, sys); % partial transpose the state
x = svd(s.data); % singular values

N = (sum(x) -1) / 2;
