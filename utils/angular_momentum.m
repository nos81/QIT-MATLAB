function J = angular_momentum(j)
% ANGULAR_MOMENTUM  Angular momentum matrices.
%  J = {Jx, Jy, Jz} = angular_momentum(j)
%
%  Returns the angular momentum matrices \vec(J)/\hbar
%  for the 2*j + 1 -dimensional subspace defined by the
%  quantum number j as a cell vector.

% Ville Bergholm 2009-2010


global qit;

n = 2*j + 1; % dimension

if (n < 1)
  error('Dimension must be one or greater.')
end

% check cache first
if (length(qit.angular_momentum) >= n && length(qit.angular_momentum{n}) > 0)
  J = qit.angular_momentum{n};
  return;
end

% raising operator in subspace J^2 = j*(j+1)
m = j;
Jplus = sparse(n,n);
for k=1:n-1
  m = m - 1;
  Jplus(k, k+1) = sqrt(j*(j+1) -m*(m+1));
end

% lowering operator
Jminus = Jplus';

% Jplus  = Jx + i*Jy
% Jminus = Jx - i*Jy

J = cell(1,3);
J{1} = 0.5*(Jplus + Jminus);
J{2} = 0.5*i*(Jminus - Jplus);
J{3} = spdiags([j:-1:-j].', 0, n, n);


% store them in the cache
qit.angular_momentum{n} = J;
