function U = single(L, t, d_in)
% SINGLE  Single-qudit operator.
%
%  U = single(L, t, d_in)
%
%  Returns the operator U corresponding to the local operator L applied
%  to subsystem t (and identity applied to the remaining subsystems).
%
%  d_in is the input dimension vector for U.

% James Whitfield 2010
% Ville Bergholm 2010-2011


if isa(L, 'lmap')
    L = L.data;
end

if d_in(t) ~= size(L, 2)
    error('Input dimensions do not match.')
end

d_out = d_in;
d_out(t) = size(L, 1);

U = lmap(op_list({{L, t}}, d_in), {d_out, d_in});
