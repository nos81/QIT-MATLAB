function varargout = heisenberg(dim, J, h, ss)
% HEISENBERG  Heisenberg model, 1D spin chain.
%  [H, dim] = heisenberg(dim, J, h)               % full Hamiltonian
%  [H, dim] = heisenberg(dim, J, h, [first last]) % partial H
%  [S_k, S_{k+1}, h_k, h_{k+1}] = heisenberg(dim, J, h, k) % operators
%
%  Returns the Hamiltonian H for an implementation of the
%  Heisenberg model, describing a one-dimensional chain of spins
%  with nearest-neighbor couplings in an external magnetic
%  field in the z direction.
%
%  The secondary calling syntax returns the coupling (S1, S2) and
%  local (h1, h2) operators just for sites k and k+1 in the chain.
%
%  dim is a vector of the dimensions of the spins, i.e. dim == [2 2 2]
%  would be a chain of three spin-1/2's.
%
%  J defines the coupling strengths in the chain. It's either a
%  3x1 vector (homogeneous couplings) or a function handle
%  J(site, component) (site-dependent couplings).
%
%  By setting J = [0 0 a] the Ising model is obtained.
%  Likewise,  J = [a a 0] gives the XY model.
%
%  h defines the couplings of the spins to the magnetic field.
%  It's either a scalar (homogeneous) or a function handle h(site)(site-dependent).
%
%  H = \sum_k (\sum_{a = x,y,z} J(k,a) Sa_{k} Sa_{k+1}) +h(k) Sz_{k})

% Ville Bergholm 2009-2010


if (nargin < 3 || nargin > 4)
  error('Usage: heisenberg(dim, J, h)')
end

if (isa(J, 'function_handle'))
  Jfunc = J;
else
  if (isvector(J) && length(J) == 3)
    Jfunc = @(k,s) J(s);
  else
    error(['J must be either a 3x1 vector or a function handle.'])
  end
end

if (isa(h, 'function_handle'))
  hfunc = h;
else
  if (isscalar(h))
    hfunc = @(k) h;
  else
    error('h must be either a scalar or a function handle.')
  end
end

if (~isvector(dim))
  error('Only 1D lattices supported for now.')
end
n = size(dim);
n = prod(n);


if (nargin == 4)
  if (isscalar(ss))
    % just return the ops for sites ss and ss+1
    S1 = angular_momentum(dim(ss));
    S2 = angular_momentum(dim(ss+1));

    for k=1:3
      varargout{1}{k} = Jfunc(ss, k) * S1{k}; % include J in S1
    end
    varargout{2} = S2;
    varargout{3} = hfunc(ss)   * S1{3}; % h1
    varargout{4} = hfunc(ss+1) * S2{3}; % h2
    return
  else
    % return partial H, ss == [start, end] (range of sites)
  end
else
  ss = [1, length(dim)]; % H for entire chain (default)
end


H = sparse(0);
S1 = angular_momentum(dim(ss(1))); % spin ops for first site

rdim = dim(ss(1):ss(2)); % reduced dimension vector
m = 1; % reduced subsystem index to go with rdim

for k=ss(1):ss(2)-1
  % coupling between spins k and k+1: S1 \cdot S2
  S2 = angular_momentum(dim(k+1)); % spin ops

  H = H +op_list({{Jfunc(k,1)*S1{1}, m; S2{1}, m+1}, {Jfunc(k,2)*S1{2}, m; S2{2}, m+1},...
                  {Jfunc(k,3)*S1{3}, m; S2{3}, m+1}, {hfunc(k)*S1{3}, m}}, rdim);

  % ops for next iteration
  S1 = S2;
  m = m+1;
end

% final term
H = H +op_list({{hfunc(ss(2)) * S1{3}, m}}, rdim);

varargout{1} = H;
varargout{2} = rdim;
