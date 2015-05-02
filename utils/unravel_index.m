function ket = unravel_index(ind, dim)
% UNRAVEL_INDEX  Convert a linear index into a multi-index array.
%  ket = unravel_index(ind, dim)
%
%  Returns the big-endian multi-index integer vector ket corresponding to
%  the linear index ind in a system with dimension vector dim.
%
%  Uses the MATLAB indexing convention (first index is 1, not zero).

% Ville Bergholm 2008-2015


d = fliplr(dim); % the big-endian convention makes this way more complicated than it should be
n = length(dim);
for k = 1:n % start from least significant digit
    %ind2sub with 2 output parms uses up the first dim given
    [ket(:,k), ind] = ind2sub(d(k:n), ind);
end
ket = fliplr(ket); % big-endian again
