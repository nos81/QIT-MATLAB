function d = dist(A, B)
% DIST  Distance between two unitary gates.
%  d = dist(A, B)
%
%  Returns \inf_{\phi \in \reals} \|A - e^{i \phi} B\|^2.

% Ville Bergholm 2007-2010


temp = A'*B;
dim = dims(temp);
d = 2*prod(dim) - 2*abs(trace(temp.data));
