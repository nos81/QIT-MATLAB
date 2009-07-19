function assert_o(a, b, tol)
% ASSERT  Octave-style assert for MATLAB.
%  assert_o(a, b, tol)
%
%  Equivalent to assert(abs(a-b) <= tol).

% Ville Bergholm 2009

assert(abs(a-b) <= tol);
