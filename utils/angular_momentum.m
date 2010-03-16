function [Jx, Jy, Jz] = angular_momentum(j)
% ANGULAR_MOMENTUM  Angular momentum matrices.
%  [Jx, Jy, Jz] = angular_momentum(j)
%
%  Returns the angular momentum matrices \vec(J)/\hbar
%  for the 2*j + 1 -dimensional subspace defined by quantum number j.

% Ville Bergholm 2009-2010


n = 2*j + 1; % dimension

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

Jx = 0.5*(Jplus + Jminus);
Jy = 0.5*i*(Jminus - Jplus);
Jz = spdiags([j:-1:-j].', 0, n, n);
