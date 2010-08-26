function U = mod_inc(x, dim, N)
% MOD_INC  Modular incrementation gate.
%  U = mod_inc(x, dim)    % N == prod(dim)
%  U = mod_inc(x, dim, N) % gate dimension prod(dim) must be >= N
%
%  Returns the gate U, which, operating on the computational state
%  |y>, increments it by x (mod N):  U |y> = |y+x (mod N)>.
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

U = sparse(dim, dim);
for y=0:N-1
  U(mod(x+y, N)+1, y+1) = 1;  % MATLAB indexing starts at 1
end

% U acts trivially for states >= N
for j=N+1:dim
  U(j,j) = 1;  % MATLAB indexing starts at 1
end
