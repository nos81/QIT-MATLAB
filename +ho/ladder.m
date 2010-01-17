function [a] = ladder(n)
% LADDER  Bosonic ladder operators.
%  a = ladder(n)
%
%  Returns the n-dimensional approximation for the bosonic
%  annihilation operator a in the number basis {|0>, |1>, ..., |n-1>}.
%  (The corresponding creation operator is of course a').

% Ville Bergholm 2009


a = diag(sqrt(1:n-1), 1);
