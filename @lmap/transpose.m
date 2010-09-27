function s = transpose(s)
% TRANSPOSE  Transpose of an lmap.
%  q = transpose(s)

% Ville Bergholm 2010


% NOTE only works as long as we have max. two indices (which we can unambiguously swap)
n = order(s);

if (n > 2)
  error('Transpose only defined for vectors and matrices.')
elseif (n == 1)
  s.dim{2} = 1; % add a singleton index
end

% swap indices
s.dim = s.dim([2 1]);

% transpose data
s.data = s.data.';
