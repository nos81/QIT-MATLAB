function [B, n_local] = tensorbasis(n, d)
% TENSORBASIS  Traceless Hermitian tensor-product basis for End(H).
%  B = tensorbasis(n, d)  % H = C_d^{\otimes n}.
%  B = tensorbasis(dim)   % H = C_{dim(1)} \otimes ... \otimes C_{dim(n)}
%
%  Returns a traceless Hermitian basis for linear operators on the Hilbert space H
%  which shares H's tensor product structure. The basis elements are tensor products
%  of Gell-Mann matrices (which in the case of qubits are equal to Pauli matrices).
%  The basis elements are normalized such that  trace(b_i'*b_j) = \delta_{ij}.
%
%  Input is either two scalars, n and d, in which case the system consists of n qu(d)its,
%  or the vector dim, which contains the dimensions of the individual subsystems.
%
%  In addition to expand Hermitian operators on H, this basis can be multiplied by
%  the imaginary unit i to be used as the antihermitian generators of SU(prod(dim)).

% Ville Bergholm 2005-2009


global qit;

if (nargin == 2)
  dim = ones(n,1) * d;
elseif (nargin == 1)
  dim = n;
  n = length(dim);
else
  error('Unknown calling syntax.')
end

% check cache first (we only cache "n qu(d)its" -type bases for convenience)
cache = false;
if (dim == dim(1))
  d = dim(1);
  if (all(size(qit.tensorbasis) >= [n, d]) && length(qit.tensorbasis{n,d}) > 0)
    B = qit.tensorbasis{n,d};
    n_local = n * (d^2 - 1);
    return;
  end
  cache = true;
end

n_elements = dim.^2;    % number of basis elements for each subsystem, incl. identity
n_all = prod(n_elements); % number of all tensor basis elements, incl. identity


local_basis{1} = [];
locals = [];
nonlocals = [];

% create the tensor basis
for k = 0:(n_all-1)  % loop over all basis elements
  temp = 1;
  sum = 0; % number of non-id. matrices included in this element
  for j=1:n  % loop over subsystems
    ind = mod(k, n_elements(j));  % which local basis element to use
    k = floor(k / n_elements(j));
    if (ind > 0)
      sum = sum + 1; % using a non-identity matrix
    end

    d = dim(j);
    if (length(local_basis) < d || isempty(local_basis{d}))
      % we need to generate the local basis for this dimension
      local_basis{d} = horzcat({eye(d)/sqrt(d)}, gellmann(d)); % hermitian, includes identity
    end

    temp = kron(temp, local_basis{d}{ind+1});  % tensor in another matrix
  end

  if (sum >= 2) % at least two non-identities => nonlocal element
    nonlocals{end+1} = temp;
  elseif (sum == 1) % skip identity
    locals{end+1} = temp;
  end
end


% arrange the elements so that local ones come first
B = horzcat(locals, nonlocals);
n_local = length(locals);
%sum(n_elements)-n  % number of local basis elements (excl. identity)


% store into cache
if (cache)
  qit.tensorbasis{n,d} = B;
end
