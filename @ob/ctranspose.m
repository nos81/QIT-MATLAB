function s = ctranspose(s)
% CTRANSPOSE  Hermitian conjugate of an lmap.

% Ville Bergholm 2010


% NOTE only works as long as we have max. two indices (which we can unambiguously swap)

% swap indices
s.dim = s.dim([2 1]);

% ctranspose data
s.data = s.data';
