function P = seq2prop(s)
% SEQ2PROP  Propagator corresponding to a control sequence.
%  P = seq2prop(s)
%
%  Returns the matrix P corresponding to the action of the control sequence s.
%  
%  Governing equation: \dot(X)(t) = (A +\sum_k u_k(t) B_k) X(t) = G(t) X(t)

%  [a, theta] ^= R_a(theta) = expm(-i*a*sigma*theta/2) = expm(-i*H*t) => H = a*sigma/2, t = theta

% Ville Bergholm 2009-2015


n = length(s.tau);
P = eye(size(s.A));

for j=1:n
  G = s.A;
  for k = 1:length(s.B)
      G = G + s.control(j, k) * s.B{k};
  end

  temp = expm(s.tau(j) * G);
  P = temp * P;
end
