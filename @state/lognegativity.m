function N = lognegativity(s, sys)
% STATE/LOGNEGATIVITY  Logarithmic negativity of the state.
%  N = lognegativity(s, sys)
%
%  Returns the logarithmic negativity of the state s wrt. the partitioning
%  given by the listing of subsystems in the vector sys.

% Ville Bergholm 2008

if (nargin < 2)
  error('Need the partitioning.')
end

N = log2(2*negativity(s, sys) + 1);
