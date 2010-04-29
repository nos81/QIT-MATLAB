function U = mod_mul(x, N, dim)
% MOD_MUL  Modular multiplication gate.
%  U = mod_mul(x, dim)    % N == dim
%  U = mod_mul(x, N, dim) % gate dimension dim must be >= N
%
%  Returns the gate U, which, operating on the computational state
%  |y>, multiplies it with x (mod N):  U |y> = |x*y (mod N)>.
%
%  If given both N and dim, U will act trivially on computational states >= N.

% Ville Bergholm 2010


if (nargin == 2)
  dim = N;
end

if (~isscalar(dim))
  dim = prod(dim); % vector of dimensions
end

if (nargin == 2)
  N = dim;
else
  if (dim < N)
    error('Register dimension must be >= N.')
  end
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
