function d = dist(A, B)
% DIST  Distance between two unitary gates.
%  d = dist(A, B)
%
%  Returns \inf_{\phi \in \reals} \|A - e^{i \phi} B\|^2.

% Ville Bergholm 2007-2009

d = 2*size(A,1) - 2*abs(trace(A'*B));
