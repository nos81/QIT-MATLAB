function U = mod_add(dim1, dim2, N)
% MOD_ADD  Modular adder gate.
%  U = mod_add(d1, d2)    % N == prod(d2)
%  U = mod_add(d1, d2, N) % target register dimension prod(d2) must be >= N
%
%  Returns the gate U, which, operating on the computational state
%  |x, y>, produces |x, y+x (mod N)>.
%  d1 and d2 are the control and target register dimensions.
%
%  If N is given, U will act trivially on target states >= N.
%
%  Notes:
%  The modular subtractor gate can be obtained by taking the
%  Hermitian conjugate of mod_add.
%  mod_add(2, 2) is equal to CNOT.

% Ville Bergholm 2010


d1 = prod(dim1);
d2 = prod(dim2);

if (nargin < 3)
  N = d2;
else
  if (d2 < N)
    error('Target register dimension must be >= N.')
  end
end

% NOTE: a real quantum computer would implement this gate using a
% sequence of reversible arithmetic gates but since we don't have
% one we might as well cheat

dim = d1 * d2;
U = sparse(dim, dim);

for a=0:d1-1
  for b=0:d2-1
    y = d2*a +b;

    if (b < N)
      x = d2*a +mod(a+b, N);
    else
      % U acts trivially for target states >= N
      x = y;
    end

    U(x+1, y+1) = 1;  % MATLAB indexing starts at 1
  end
end

dim = [dim1, dim2];
U = lmap(U, {dim, dim});
