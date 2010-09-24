function s = transpose(s)
% TRANSPOSE  Transpose of an lmap.

% Ville Bergholm 2010


% NOTE only works as long as we have max. two indices (which we can unambiguously swap)

% swap indices
s.dim = s.dim([2 1]);

% transpose data
s.data = s.data.';
