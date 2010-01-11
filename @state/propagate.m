function out = propagate(s, H, t, varargin)
% STATE/PROPAGATE  Propagate the state in time using a Hamiltonian or a Liouvillian.
%
%  out = propagate(s, H, t [, out_func, odeopts])
%  out = propagate(s, L, t [, out_func, odeopts])
%
%  Propagates the state s using the Hamiltonian H for the time t,
%  returns the resulting state.
%
%  H can either be a matrix or, for time-dependent Hamiltonians, a function handle.
%  The H function must take a time instance as input and return the corresponding H matrix.
%
%  If t is a vector of increasing time instances, returns a cell array
%  containing the propagated state for each time given in t.
%
%  Optional parameters (can be given in any order):
%    out_func: If given, for each time instance propagate returns out_func(s(t), H(t)).
%    odeopts:  Options struct for MATLAB ODE solvers from the odeset function.
%
%  out == expm(-i*H*t)*|s>
%  out == inv_vec(expm(L*t)*vec(\rho_s))


% Ville Bergholm 2008-2010
% James Whitfield 2009


if (nargin < 3)
  error('Needs a state, a Hamiltonian and a time.');
end

out_func    = @(x,h) x; % if no out_func is given, use a NOP
odeopts = [];
t_dependent = false; % time independent Hamiltonian by default

if (isnumeric(H))
  % time independent
  d = size(H, 2);
elseif (isa(H, 'function_handle'))
  % time dependent
  t_dependent = true;
  d = size(H(0), 2);
else
  error('The Hamiltonian parameter has to be either a matrix or a function handle.')
end

% process optional arguments
for k=1:nargin-3
  switch class(varargin{k})
    case 'function_handle'
      out_func = varargin{k};

    case 'struct'
      odeopts = varargin{k};

    otherwise
      error('Unknown optional parameter type.');
  end
end


n = length(t); % number of time instances we are interested in
out = cell(1, n);

temp = size(s.data, 1); % system dimension

if (d == temp)
  % Hamiltonian evolution
  superop = false;
elseif (d == temp^2)
  % Liouvillian superoperator evolution
  superop = true;

  s = to_op(s);
  s0 = vec(s.data); % state operator, arranged as a col vector
else
  error('Dimension of the Hamiltonian does not match the dimension of the state.');
end


% mixed states require a little extra effort
ddd = size(s.data);
function d = mixed_fun(t, y, H)
  y = reshape(y, ddd); % into a matrix
  d = -i * (H(t) * y  -y * H(t)); % derivative for the solver
  d = d(:); % back into a vector
end

if (t_dependent)
  % time dependent case, use ODE solver

  if (superop == false)
    s0 = s.data;
    if (ddd(2) == 1)
      % pure
      odefun = @(t, y, H) -i * H(t) * y; % derivative function for the solver
    else
      % mixed
      odefun = @mixed_fun;
    end
  else
    odefun = @(t, y, H) H(t) * y;
  end

  skip = 0;
  if (t(1) ~= 0)
    t = [0, t]; % ODE solver needs to be told that t0 = 0
    skip = 1;
  end

  if (length(t) < 3)
    t = [t, t(end)+1e5*eps]; % add another time point to get reasonable output from solver
  end

  %odeopts = odeset(odeopts, 'OutputFcn', @(t,y,flag) odeout(t, y, flag, H));

  [t_out, s_out] = ode45(odefun, t, s0, odeopts, H);
  s_out = s_out.'; % convention: states in columns

  % apply out_func
  for k=1:n
    if (superop == false)
      s.data = reshape(s_out(:, k+skip), ddd);
    else
      s.data = inv_vec(s_out(:, k+skip));
    end
    out{k} = out_func(s, H(t_out(k+skip)));
  end

else
  % time independent case

  if (superop == false)
    % eigendecomposition
    [v, d] = eig(H); 
    d = diag(d);
    for k=1:n
      U = v * diag(exp(-i * d * t(k))) * v';
      out{k} = out_func(u_propagate(s, U), H);
      %out{k} = out_func(u_propagate(s, expm(-i*H*t(k))), H);
    end

  else
    % Krylov subspace method
    [w, err] = expv(t, H, s0);
    for k=1:n
      s.data = inv_vec(w(:,k));
      out{k} = out_func(s, H);
    end
  end
end

if (n == 1)
  % single output, don't bother with a list
  out = out{1};
end

end


%function status = odeout(t, y, flag, H)
%if isempty(flag)
%  sdfsd
%end
%status = 0;
%end


function obsolete_stuff()
  T = t(end); % final time
  dt = T/steps; % maximum step size
  r  = 1; % index into t
  tt = 0; % time now
  h = H(0); % value of Liouvillian now

  while (tt < T)
    step = min(dt, t(r)-tt);

    % propagate, advance time
    [svec, err] = expv(step, h, svec);
    s = u_propagate(s, expm(-i*h*step));

    tt = tt + step;
    h = H(tt); % value of Hamiltonian now

    % should we output something?
    if (t(r) <= tt)
      s.data = inv_vec(svec);
      out{r} = out_func(s, h);
      r = r + 1; % next output
    end
  end
end
