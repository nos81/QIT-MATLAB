function U = single(L, targ, dim)
% SINGLE  Single-qudit operator.
%
%  U = single(L, targ, dim)
%
%  Returns the operator U corresponding to the local operator L applied
%  to subsystem targ (and identity applied to the remaining subsystems).
%
%  dim is the dimension vector for U.

% James Whitfield 2010
% Ville Bergholm 2010


d1 = dim;
d1(targ) = prod(size(L.data, 1));
U = lmap(op_list({{L.data, targ}}, dim), {d1, dim});
