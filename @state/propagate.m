function out = propagate(s, H, t, varargin)
% STATE/PROPAGATE  Propagate a quantum state in time using a Hamiltonian or a Liouvillian.
%
%  out = propagate(s, H, t [, out_func])
%  out = propagate(s, H_func, t, steps [, out_func])
%
%  Propagates the state s using the Hamiltonian H for the time t,
%  returns the resulting state.
%  If t is a vector, returns a cell array containing the propagated
%  state for each time point given in t.
%
%  out == expm(-i*H*t)*|s>
%  out == inv_vec(expm(L*t)*vec(\rho_s))


% Ville Bergholm 2008-2009
% James Whitfield 2009


if (nargin < 3)
  error('Needs a state, Hamiltonian and time duration(s).');
end

out_func    = @(x,h)x; % if no out_func is given, use a NOP
t_dependent = false; % time independent Hamiltonian by default

if (isnumeric(H))
  % time independent, no steps parameter
  d = size(H, 2);
  if (nargin == 4)
    out_func = varargin{1};
  end

elseif (isa(H, 'function_handle'))
  % time dependent
  t_dependent = true;
  d = size(H(0), 2);

  if (nargin < 4)
    error('Number of steps needed.');
  else
    steps = varargin{1};
    if (nargin == 5)
      out_func = varargin{2};
    end
  end

else
  error('The Hamiltonian parameter has to be either a matrix or a function handle.')
end


n = length(t); % number of time instances we are interested in
temp = size(s.data, 1); % system dimension

if (d == temp) % Hamiltonian evolution
  if (~t_dependent)
    if (n == 1)
      % single output, don't bother with a list
      out = out_func(u_propagate(s, expm(-i*H*t)), H);
      return
    end

    % several time instances requested, use eigendecomposition
    out = cell(n, 1);

    [v, d] = eig(H); 
    d = diag(d);
    for k=1:n
      U = v*diag(exp(-i*d*t(k)))*v';
      out{k} = out_func(u_propagate(s, U), H);
    end

  else
    T = t(end); % final time
    dt = T/steps; % maximum step size
    r  = 1; % index of next t we are interested in
    tt = 0; % time now
    h = H(0); % value of Hamiltonian now

    while (tt < T)
      step = min(dt, t(r)-tt);

      % propagate, advance time
      s = u_propagate(s, expm(-i*h*step));
      tt = tt + step;
      h = H(tt); % value of Hamiltonian now

      % should we output something?
      if (t(r) <= tt)
	out{r} = out_func(s, h);
	r = r + 1; % next output
      end
    end
  end

elseif (d == temp^2) % Liouvillian superoperator evolution

  svec = vec(to_op(s).data); % state operator, arranged as a col vector

  if (~t_dependent)
    if (n == 1)
      % single output, don't bother with a list
      [svec, err] = expv(t, H, svec);
      s.data = inv_vec(svec);
      out = out_func(s, H);
      return
    end

    % several time instances requested
    out = cell(n, 1);

    [w, err] = expv(t, H, svec);
    for k=1:n
      s.data = inv_vec(w(:,k));
      out{k} = out_func(s, H);
    end

  else
    T = t(end); % final time
    dt = T/steps; % maximum step size
    r  = 1; % index into t
    tt = 0; % time now
    h = H(0); % value of Liouvillian now

    while (tt < T)
      step = min(dt, t(r)-tt);

      % propagate, advance time
      [svec, err] = expv(step, h, svec);

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

else
  error('Dimension of the Hamiltonian does not match the dimension of the state.');
end
