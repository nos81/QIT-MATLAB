function out = propagate(s, H, t, varargin)
% PROPAGATE  Propagate the state continuously in time.
%
%  out = propagate(s, H, t [, out_func, odeopts])
%  out = propagate(s, L, t [, out_func, odeopts])
%  out = propagate(s, {H, {A_i}}, t [, out_func, odeopts])
%
%  Propagates the state s using the given generator(s) for the time t,
%  returns the resulting state.
%
%  The generator can either be a Hamiltonian matrix H or, for time-dependent
%  Hamiltonians, a function handle H(t) which takes a time instance t
%  as input and return the corresponding H matrix.
%
%  Alternatively, the generator can also be a Liouvillian superoperator, or
%  a list consisting of a Hamiltonian and a list of Lindblad operators.
%
%  If t is a vector of increasing time instances, returns a cell array
%  containing the propagated state for each time given in t.
%
%  Optional parameters (can be given in any order):
%    out_func: If given, for each time instance propagate returns out_func(s(t), H(t)).
%    odeopts:  Options struct for MATLAB ODE solvers from the odeset function.

%  out == expm(-i*H*t)*|s>
%  out == inv_vec(expm(L*t)*vec(\rho_s))

% Ville Bergholm 2008-2010
% James Whitfield 2009


if (nargin < 3)
  error('Needs a state, a generator and a time.');
end

out_func = @(x,h) x; % if no out_func is given, use a NOP

odeopts = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'Vectorized', 'on');

n = length(t); % number of time instances we are interested in
out = cell(1, n);
dim = size(s.data); % system dimension

if (isa(H, 'function_handle'))
  % time dependent
  t_dependent = true;
  F = H;
  H = F(0);
else
  % time independent
  t_dependent = false;
end

if (isnumeric(H))
  % matrix
  dim_H = size(H, 2);

  if (dim_H == dim(1))
    gen = 'H';
  elseif (dim_H == dim(1)^2)
    gen = 'L';
    s = to_op(s);
  else
    error('Dimension of the generator does not match the dimension of the state.');
  end
  
elseif (iscell(H))
  % list of Lindblad operators
  dim_H = size(H{1}, 2);
  if (dim_H == dim(1))
    gen = 'A';
    s = to_op(s);

    % HACK, in this case we use ode45 anyway
    if (~t_dependent)
      t_dependent = true; 
      F = @(t) H; % ops stay constant
    end
  else
    error('Dimension of the Lindblad ops does not match the dimension of the state.');
  end

else
  error(['The second parameter has to be either a matrix, a cell array, '...
         'or a function handle that returns a matrix or a cell array.'])
end

dim = size(s.data); % may have been switched to operator representation


% process optional arguments
for k=1:nargin-3
  switch class(varargin{k})
    case 'function_handle'
      out_func = varargin{k};

    case 'struct'
      odeopts = odeset(odeopts, varargin{k});

    otherwise
      error('Unknown optional parameter type.');
  end
end


% derivative functions for the solver

function d = lindblad_fun(t, y, F)
  X = F(t);
  A = X{2};
  A = A(:);

  d = zeros(size(y));
  % lame vectorization
  for loc1_k=1:size(y, 2)
    x = reshape(y(:,loc1_k), dim); % into a matrix
  
    % Hamiltonian
    temp = -1i * (X{1} * x  -x * X{1});
    % Lindblad operators
    for j=1:length(A)
      ac = A{j}'*A{j};
      temp = temp +A{j}*x*A{j}' -0.5*(ac*x + x*ac);
    end
    d(:,loc1_k) = temp(:); % back into a vector
  end
end

function d = mixed_fun(t, y, F)
  H = F(t);
  
  d = zeros(size(y));
  % vectorization
  for loc2_k=1:size(y, 2)
    x = reshape(y(:,loc2_k), dim); % into a matrix
    temp = -1i * (H * x  -x * H);
    d(:,loc2_k) = temp(:); % back into a vector
  end
end


if (t_dependent)
  % time dependent case, use ODE solver

  switch (gen)
    case 'H'
      % Hamiltonian
      if (dim(2) == 1)
        % pure state
        odefun = @(t, y, F) -1i * F(t) * y; % derivative function for the solver
      else
        % mixed state
        odefun = @mixed_fun;
      end

    case 'L'
      % Liouvillian
      odefun = @(t, y, F) F(t) * y;
      %odeopts = odeset(odeopts, 'Jacobian', F);

    case 'A'
      % Hamiltonian and Lindblad operators in a list
      odefun = @lindblad_fun;
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

  [t_out, s_out] = ode45(odefun, t, s.data, odeopts, F);
  % s_out has states in columns, row i corresponds to t(i)

  % apply out_func
  for k=1:n
    % this works because ode45 automatically expands input data into a col vector
    s.data = inv_vec(s_out(k+skip,:), dim);
    out{k} = out_func(s, F(t_out(k+skip)));
  end

else
  % time independent case

  switch (gen)
    case 'H'
      if (length(H) < 500)
        % eigendecomposition
        [v, d] = eig(full(H)); % TODO eigs?
        d = diag(d);
        for k=1:n
          U = v * diag(exp(-1i * t(k) * d)) * v';
          out{k} = out_func(s.prop(U), H);
          %out{k} = out_func(s.prop(expm(-i*H*t(k))), H);
        end
      else
        % Krylov subspace method
        [w, err] = expv(-1i*t, H, s.data);
        for k=1:n
          s.data = w(:,k);
          out{k} = out_func(s, H);
        end
      end
      
    case 'L'
      % Krylov subspace method
      [w, err] = expv(t, H, vec(s.data));
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
