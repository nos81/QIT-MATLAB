function varargout = spin1d(dim, J, h, ss)
% SPIN1D  Spin model, 1D spin chain.
%  [H, dim] = spin1d(dim, J, h)               % full Hamiltonian
%  [H, dim] = spin1d(dim, J, h, [first last]) % partial H
%  [S_k, S_{k+1}, h_k, h_{k+1}] = spin1d(dim, J, h, k) % operators
%
%  Returns the Hamiltonian H for an implementation of the
%  Hamiltonian model, describing a one-dimensional chain of spins
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
%  upper triangular 3x3 matrix (homogeneous couplings) or a function handle
%  J(site) that returns upper triangular matrices (site-dependent couplings).
%
%  By setting J = [0 0 0;
% 		   0 0 0;
% 		   0 0 a] the Ising model is obtained.
%  Likewise, J =  [a 0 0;
% 		   0 a 0;
% 		   0 0 0] yields the XY model
%
%  h defines the couplings of the spins to the magnetic field and is a 3x1 vector
%  (homogeneous) or a function handle h(site) that returns 
%  3x1 vectors (site-dependent).
%
%  H = \sum_k (\sum{a=x,y,z}\sum{b=x,y,z} J_k(a,b) S_k^a S_k^b + h_k(a)S_k^a)

% Ville Bergholm 2009-2010
% James Whitfield 2010


if (nargin < 3 || nargin > 4)
  error('Usage: spin1d(dim, J, h)')
end

if (isa(J, 'function_handle'))
  Jfunc = J;
else
  if (size(J,1)==size(J,2) && length(J) == 3)
    Jfunc = @(k) J;
  else
    error(['J must be either a 3x3 matrix or a function handle.'])
  end
end

if (isa(h, 'function_handle'))
  hfunc = h;
else
  if (isvector(h))
    hfunc = @(k) h;
  else
    error('h must be either a vector or a function handle.')
  end
end

if (~isvector(dim))
  error('Only 1D spin chains.')
end
n = size(dim);
n = prod(n);


if (nargin == 4)
  if (isscalar(ss))
    % just return the ops for sites ss and ss+1
    S1 = angular_momentum(dim(ss));
    S2 = angular_momentum(dim(ss+1));
    J=Jfunc(ss);
    %matJ=
	%[XX,XY,XZ ]
	%[   YY,YZ ]
	%[      ZZ ]
    h=hfunc(ss);
    %h=
	%[X]
	%[Y]
	%[Z]
    idx=1;
    for comp1=1:3
   	 for comp2=(1+(comp1-1)):3
		if J(comp1,comp2)
			varargout{1}{idx} = J(comp1,comp2) * S1{comp1}; % include J in S1
			varargout{2}{idx} = S2{comp2};
			idx=idx+1;
		end
    	end

	varargout{3}{comp1} = h(comp1) * S1{comp1}; % h1
	varargout{4}{comp1} = h(comp1) * S2{comp1}; % h1
    end

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
  J=Jfunc(k);
  h=hfunc(k);
  %H = H +op_list({{Jfunc(k,1)*S1{1}, m; S2{1}, m+1}, {Jfunc(k,2)*S1{2}, m; S2{2}, m+1},...
  %                {Jfunc(k,3)*S1{3}, m; S2{3}, m+1}, {hfunc(k)*S1{3}, m}}, rdim);

  H = H+op_list({...
  {J(1,1)*S1{1},m; S2{1},m+1},{J(1,2)*S1{1}, m; S2{2}, m+1},{J(1,3)*S1{1}, m; S2{3}, m+1},...
                              {J(2,2)*S1{2}, m; S2{2}, m+1},{J(2,3)*S1{2}, m; S2{3}, m+1},...
                                                            {J(3,3)*S1{3}, m; S2{3}, m+1},...
  {h(1)*S1{1}, m},            {h(2)*S1{2}, m},              {h(3)*S1{3},m}},rdim);

  % ops for next iteration
  S1 = S2;
  m = m+1;
end

% final term
h=hfunc(ss(2));
H = H +op_list({{h(1) * S1{1}, m},{h(2) * S1{2}, m},{h(3) * S1{3}, m}}, rdim);

varargout{1} = H;
varargout{2} = rdim;
