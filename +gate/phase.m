function U = phase(theta)
% PHASE  Phase shift gate.
%  U = phase(theta)
%
%  Returns the phase shift gate U = diag([1, e^(i*theta)]).

% Ville Bergholm 2010


n = length(theta)+1; % first phase is implicitly taken to be 1
U = spdiags([1, exp(i*theta)].', 0, n, n);
