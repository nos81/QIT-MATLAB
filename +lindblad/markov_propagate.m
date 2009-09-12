function out = markov_propagate(s, H, D, b, t, func)
% LINDBLAD/MARKOV_PROPAGATE  Propagate a system in time using a Hamiltonian and a coupled bath.
%
%  out = markov_propagate(s, H, D, b, t [, func])
%
%  Propagates the state s using the Hamiltonian H for the time t,
%  under the influence of the bath b coupled to s through the operator D,
%  returns the resulting state.
%  If t is a vector, returns a cell array containing the propagated
%  state for each time point given in t.
%
%  If a function handle func is given, each propagated state is inputted
%  to it and the result is given as output.

% Ville Bergholm 2009


if (nargin < 6)
    func = @(x) x; % no function given, use a NOP

    if (nargin < 5)
      error('Needs a state, a Hamiltonian, bath coupling, bath, and a time duration.');
    end
end

if (size(s.data, 1) ~= size(H, 2))
  error('Dimension of the Hamiltonian does not match the dimension of the state.');
end


n = length(t);
out = cell(n, 1);
    
L = liouvillian(H, D, b); % build noisy liouvillian

s = operator(s);
svec = vec(s.data); % must be state op

switch 'expv'
case 'expm'
  for k=1:n
      s.data = inv_vec(expm(L * t(k))*svec);
      out{k} = func(s);
  end

case 'eig'
  % slow, inaccurate in octave
  [v, d] = eig(L); 
  d = diag(d);
  %norm(v*v'-eye(256))
  plot(d, 'rx')
  for k=1:n
      U = v*diag(exp(d * t(k)))*v';
      s.data = inv_vec(U*svec);
      out{k} = func(s);
  end

case 'expv'
  % Krylov subspace method from Expokit
  [w, err] = expv(t, L, svec);
  for k=1:n
      s.data = inv_vec(w(:,k));
      out{k} = func(s);
  end
end
end
