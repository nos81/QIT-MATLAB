function out = propagate(s, H, t, func)
% STATE/PROPAGATE  Propagate a quantum state in time using a Hamiltonian.
%  out = propagate(s, H, t [, func])
%
%  Propagates the state s using the Hamiltonian H for the time t,
%  returns the resulting state.
%  If t is a vector, returns a cell array containing the propagated
%  state for each time point given in t.
%
%  If a function handle func is given, each propagated state is inputted
%  to it and the result is given as output.

% Ville Bergholm 2008-2009
% James Whitfield 2009


if (nargin < 4)
    func = @(x) x; % no function given, use a NOP

    if (nargin < 3)
      error('Needs a state, a Hamiltonian and a time duration.');
    end
end

if (size(s.data, 1) ~= size(H, 2))
  error('Dimension of the Hamiltonian does not match the dimension of the state.');
end

n = length(t);

if (n == 1)
  out = func(u_propagate(s, expm(-i*H*t)));
  return
end

% several time instances requested, use eigendecomposition
out = cell(n, 1);
[v, d] = eig(H); 

for k=1:n
  U = v*diag(exp(-i*diag(d)*t(k)))*v';
  out{k} = func(u_propagate(s, U));
end
