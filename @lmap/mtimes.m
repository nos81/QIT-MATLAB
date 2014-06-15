function b = mtimes(a, b)
% MTIMES  Multiplication of lmaps by lmaps and scalars.
%  s = mtimes(a, b)
%
%  Returns the product of a and b.

% Ville Bergholm 2010-2014


if isa(a, 'lmap') && isa(b, 'lmap')
  % in lieu of a full contraction method, we only multiply vector and matrices

  is_concatenable(a, b);
  b = lmap(a.data * b.data, {a.dim{1}, b.dim{2}});

elseif isscalar(a) && isnumeric(a)
  b.data = a * b.data;
elseif isscalar(b) && isnumeric(b)
  temp = b;
  b = a;
  b.data = temp * b.data;
else
  error('lmaps can only be multiplied by scalars and lmaps.')
end
