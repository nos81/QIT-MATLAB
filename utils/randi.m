function a = randi(n)
% RANDI  Uniformly distributed pseudorandom integers.
%  a = randi(n)
%
% Returns a pseudorandom integer uniformly distributed between 1 and n.

% Ville Bergholm 2010


  a = ceil(rand*n);
end
