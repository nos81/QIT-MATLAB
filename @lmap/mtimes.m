function b = mtimes(a, b)
% MTIMES  Multiplication of lmaps by lmaps and scalars.
%  s = times(a, b)
%
%  Returns product of a and b.

% Ville Bergholm 2010


if (isa(a, 'lmap') && isa(b, 'lmap'))
  if (~isequal(a.dim{2}, b.dim{1}))
    error('The dimensions do not match.')
  end

  b.data = a.data*b.data;
  b.dim{1} = a.dim{1};

elseif (isscalar(a) && isnumeric(a))
  b.data = a * b.data;
elseif (isscalar(b) && isnumeric(b))
  temp = b;
  b = a;
  b.data = temp * b.data;
else
  error('Lmaps can only be multiplied by scalars and lmaps.')
end
