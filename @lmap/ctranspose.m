function s = ctranspose(s)
% CTRANSPOSE  Hermitian conjugate of an lmap.
%  q = ctranspose(s)

% Ville Bergholm 2010


% NOTE only works as long as we have max. two indices (which we can unambiguously swap)
n = order(s);

if (n > 2)
  error('Hermitian conjugate only defined for vectors and matrices.')
elseif (n == 1)
  s.dim{2} = 1; % add a singleton index
end

% swap indices
s.dim = s.dim([2 1]);

% ctranspose data
s.data = s.data';
