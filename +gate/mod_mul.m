function U = mod_mul(x, dim, N)
% MOD_MUL  Modular multiplication gate.
%  U = mod_mul(x, dim)    % N == prod(dim)
%  U = mod_mul(x, dim, N) % gate dimension prod(dim) must be >= N
%
%  Returns the gate U, which, operating on the computational state
%  |y>, multiplies it by x (mod N):  U |y> = |x*y (mod N)>.
%  x and N must be coprime for the operation to be reversible.
%
%  If N is given, U will act trivially on computational states >= N.

% Ville Bergholm 2010


if (~isscalar(dim))
  dim = prod(dim); % vector of dimensions
end

if (nargin < 3)
  N = dim;
else
  if (dim < N)
    error('Gate dimension must be >= N.')
  end
end

if (gcd(x, N) ~= 1)
  error('x and N must be coprime for the mul operation to be reversible.')
end

% NOTE: a real quantum computer would implement this gate using a
% sequence of reversible arithmetic gates but since we don't have
% one we might as well cheat
U = sparse(dim, dim);
for j=1:N
  y = j-1;
  temp = mod(x*y, N);
  U(temp+1,j) = 1;
end

% U acts trivially for states >= N
for j=N+1:dim
  U(j,j) = 1;
end
