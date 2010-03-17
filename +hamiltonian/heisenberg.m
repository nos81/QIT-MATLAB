function [H, dim, varargout] = heisenberg(dim, J, h, site)
% HEISENBERG  Heisenberg model, 1D spin chain.
%  [H, dim] = heisenberg(dim, J, h) % full Hamiltonian
%  [S_k, S_{k+1}, h_k, h_{k+1}] = heisenberg(dim, J, h, k) % individual ops
%
%  Returns the Hamiltonian H for an implementation of the
%  Heisenberg model, describing a one-dimensional chain of spins
%  with nearest-neighbor couplings in an external magnetic
%  field in the z direction.
%
%  The second calling syntax returns the coupling and local
%  operators just for sites k and k+1 in the chain.
%
%  dim is a vector of the dimensions of the spins, i.e. [2 2 2]
%  would be a chain of three spin-1/2's.
%
%  J defines the coupling strengths in the chain. It's either a
%  3x1 vector (homogeneous couplings) or a function handle
%  (site-dependent couplings).
%
%  By setting J = [0 0 a] the Ising model is obtained.
%  Likewise,  J = [a a 0] gives the XY model.
%
%  h gives the couplings of the spins to the magnetic field.
%  It's either a scalar (homogeneous) or a function handle (site-dependent).
%
%  H = \sum_k (\sum_{a = x,y,z} J(k,a) Sa_{k} Sa_{k+1}) +h_{k} Sz_{k})

% Ville Bergholm 2009-2010


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


% spin ops
S1 = angular_momentum((dim(1)-1)/2); % from dimension to spin quantum number

if (nargin == 4)
  % just return the ops for site and site+1
  % dim has the corresponding two dimensions
  S2 = angular_momentum((dim(2)-1)/2);
  for k=1:3
    H{k} = Jfunc(site, k) * S1{k}; % include J in S1
  end

  % FIXME ugly hack
  H = S1;
  dim = S2;
  varargout{1} = hfunc(site)   * S1{3}; % h1
  varargout{2} = hfunc(site+1) * S2{3}; % h2
  return
end


H = sparse(0);

for k=1:n-1
  % coupling between spins k and k+1: S1 \cdot S2
  S2 = angular_momentum((dim(k+1)-1)/2); % spin ops

  H = H +op_list({{Jfunc(k,1)*S1{1}, k; S2{1}, k+1}, {Jfunc(k,2)*S1{2}, k; S2{2}, k+1},...
                  {Jfunc(k,3)*S1{3}, k; S2{3}, k+1}, {hfunc(k)*S1{3}, k}}, dim);

  % ops for next iteration
  S1 = S2;
end

% final term
H = H +op_list({{hfunc(n) * S1{3}, n}}, dim);
