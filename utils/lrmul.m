function S = lrmul(L, R)
% LRMUL  Superoperator equivalent for multiplying both from left and right.
%  S = lrmul(L, R);
%
%  L*rho*R == inv_vec(lrmul(L, R)*vec(rho))

% Ville Bergholm 2009


if (issquare(L) && issquare(R))
  S = kron(R.', L); % simplifies to this when L and R are both square
end

n = size(L, 1);
q = size(R, 1);

S = kron(R.', eye(n)) * kron(eye(q), L);
end

function b = issquare(A)
  s = size(A);
  b = (s(1) == s(2));
end