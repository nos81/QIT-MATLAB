function x = mpower(a, b)
% MPOWER  Exponentiation of lmaps by scalars.
%  x = mpower(a, b)
%
%  Returns a^b.

% Ville Bergholm 2010


if (isa(a, 'lmap') && isscalar(b) && isnumeric(b))

  if (floor(b) ~= b)
    error('Can only handle positive integer exponents for now.')
  end

  n = order(a);
  if (n ~= 2)
    error('a is not a matrix.');
  end

  if (~isequal(a.dim{2}, a.dim{1}))
    error('The dimensions do not match.')
  end
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
