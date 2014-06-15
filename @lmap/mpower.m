function x = mpower(a, b)
% MPOWER  Exponentiation of lmaps by scalars.
%  x = mpower(a, b)
%
%  Returns a^b.

% Ville Bergholm 2010-2014


if isa(a, 'lmap') && isscalar(b) && isnumeric(b)

  if floor(b) ~= b
    error('Can only handle positive integer exponents for now.')
  end

  is_concatenable(a, a);
  dim = a.dim{1};

  % exponentiation by repeated squaring
  b = fliplr(dec2bin(b) - '0'); % exponent in little-endian binary
  x = lmap(speye(prod(dim)), {dim, dim});
  for k = 1:length(b)
    if (b(k))
      x = a * x;
    end
    a = a*a; % square it
  end

else
  error('lmaps can only be exponentiated by scalars.')
end
