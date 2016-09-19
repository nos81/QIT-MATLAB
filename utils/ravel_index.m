function ind = ravel_index(ket, dim)
% RAVEL_INDEX  Convert a multi-index array into a linear index.
%  ind = ravel_index(ket, dim)
%
%  Returns the linear index ind corresponding to the big-endian
%  multi-index integer vector ket in a system with dimension vector dim.
%
%  Uses the MATLAB indexing convention (the first index value is 1, not zero).

% Ville Bergholm 2016

n = length(ket);
ind = 0;

if any(ket > dim) || any(ket < 1)
    error('Index value out of range.')
end

d = fliplr(cumprod(fliplr(dim)));
d = [d(2:end), 1];
ind = dot(ket-1, d)+1;
