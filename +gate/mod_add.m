function U = mod_add(x, N, dim)
% MOD_ADD  Modular addition gate.
%  U = mod_add(x, dim)    % N == dim
%  U = mod_add(x, N, dim) % gate dimension dim must be >= N
%
%  Returns the gate U, which, operating on the computational state
%  |y>, increments it by x (mod N):  U |y> = |y+x (mod N)>.
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
for y=1:N
  U(mod(x+y-1, N)+1, y) = 1; % MATLAB indexing starts at 1
end

% U acts trivially for states >= N
for j=N+1:dim
  U(j,j) = 1;
end
