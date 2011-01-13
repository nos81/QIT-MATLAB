function b = mtimes(a, b)
% MTIMES  Multiplication of lmaps by lmaps and scalars.
%  s = mtimes(a, b)
%
%  Returns product of a and b.

% Ville Bergholm 2010


if (isa(a, 'lmap') && isa(b, 'lmap'))
  % in lieu of a full contraction method, we only multiply vector and matrices

  n = order(a);
  m = order(b);

  if (n ~= 2)
    error('a is not a matrix.');
  end
  if (m ~= 2)
    error('b is not a matrix.');
  end

  if (~isequal(a.dim{2}, b.dim{1}))
    error('The dimensions do not match.')
  end

  b = lmap(a.data*b.data, {a.dim{1}, b.dim{2}});

elseif (isscalar(a) && isnumeric(a))
  b.data = a * b.data;
elseif (isscalar(b) && isnumeric(b))
  temp = b;
  b = a;
  b.data = temp * b.data;
else
  error('lmaps can only be multiplied by scalars and lmaps.')
end
