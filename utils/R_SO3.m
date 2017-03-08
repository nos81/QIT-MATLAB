function R = R_SO3(a)
% R_SO3  SO(3) rotation around a given vector.
%  O = R_SO3(a)
%
%  Returns the SO(3) rotation about the unit three-vector a/|a| by the angle |a|.

% Ville Bergholm 2016

ox = [0, 0, 0; 0, 0, -1; 0, 1, 0];
oy = [0, 0, 1; 0, 0, 0; -1, 0, 0];
oz = [0, -1, 0; 1, 0, 0; 0, 0, 0];

R = expm(a(1)*ox +a(2)*oy +a(3)*oz);
