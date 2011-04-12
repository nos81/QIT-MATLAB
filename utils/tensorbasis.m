function [B, local] = tensorbasis(n, d)
% TENSORBASIS  Hermitian tensor-product basis for End(H).
%  B = tensorbasis(n, d)  % H = C_d^{\otimes n}.
%  B = tensorbasis(dim)   % H = C_{dim(1)} \otimes ... \otimes C_{dim(n)}
%
%  Returns a Hermitian basis for linear operators on the Hilbert space H
%  which shares H's tensor product structure. The basis elements are tensor products
%  of Gell-Mann matrices (which in the case of qubits are equal to Pauli matrices).
%  The basis elements are normalized such that \trace(b_i' * b_j) = \delta_{ij}.
%
%  Input is either two scalars, n and d, in which case the system consists of n qu(d)its,
%  or the vector dim, which contains the dimensions of the individual subsystems.
%
%  In addition to expanding Hermitian operators on H, this basis can be multiplied by
%  the imaginary unit i to obtain the antihermitian generators of U(prod(dim)).

% Ville Bergholm 2005-2011


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
  if (all(size(qit.tensorbasis) >= [n, d]) && ~isempty(qit.tensorbasis{n,d}))
    B = qit.tensorbasis{n,d};
    local = qit.tensorbasislocal{n,d};
    return;
  end
  cache = true;
end

n_elements = dim.^2;      % number of basis elements for each subsystem, incl. identity
n_all = prod(n_elements); % number of all tensor basis elements, incl. identity

B = cell(1, n_all); % basis
local = false(1, n_all); % logical array, is the corresponding basis element local?

% create the tensor basis
for k = 0:(n_all-1)  % loop over all basis elements
  temp = 1; % basis element being built
  rem = k; % remainder
  sum = 0; % number of non-id. matrices included in this element
  for j=1:n  % loop over subsystems
    ind = mod(rem, n_elements(j));  % which local basis element to use
    rem = floor(rem / n_elements(j));
    if (ind > 0)
      sum = sum + 1; % using a non-identity matrix
    end

    d = dim(j);
    L = horzcat({eye(d)/sqrt(d)}, gellmann(d)); % Gell-Mann basis for the subsystem (cached)
    temp = kron(temp, L{ind+1});  % tensor in another matrix
  end

  B{k+1} = temp;  
  local(k+1) = (sum < 2); % at least two non-identities => nonlocal element
end


% store into cache
if (cache)
  qit.tensorbasis{n,d} = B;
  qit.tensorbasislocal{n,d} = local;
end
