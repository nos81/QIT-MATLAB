function U = seq2prop(s)
% SEQ2PROP  SU(2) propagator corresponding to a single-qubit control sequence.
%  U = seq2prop(s)
%
%  Returns the SU(2) rotation matrix U corresponding to the
%  action of the single-qubit control sequence s alone.

%  [a, theta] ^= R_a(theta) = expm(-i*a*sigma*theta/2) = expm(-i*H*t) => H = a*sigma/2, t = theta

% Ville Bergholm 2009


global qit;

n = size(s, 1);
U = eye(2);

for k=1:n
  H = 0.5*(qit.sx*s(k,1) +qit.sy*s(k,2) +qit.sz*s(k,3));
  t = s(k,end);
  temp = expm(-i*H*t);
  U = temp * U;
end
