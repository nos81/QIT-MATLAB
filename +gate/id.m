function U = id(dim)
% ID  Identity gate.
%  I = id(dim)
%
%  Returns the identity gate I for the specified system.
%  dim is a vector of subsystem dimensions.

% Ville Bergholm 2010


U = lmap(speye(prod(dim)), {dim, dim});
