function s = mps(A)
% MPS  Matrix product states.
%  s = mps(A)
%
%  Returns the matrix product state s corresponding to the
%  projection matrices in the cell vector A.
%
%  The length of A determines the number of subsystems in s.
%  Each element A{k} is a cell vector of matrices, the number of
%  which determines the dimension of subsystem k.
%
%  Using the syntax in the reference A{k}{i} == A^{[k]}_i.

%! D. Perez-Garcia et al., "Matrix Product State Representations", Quant. Inf. Comput. 7, 401 (2007). arXiv.org:quant-ph/0608197
% Ville Bergholm 2010


if (~iscell(A))
  error('A must be a cell vector.')
end

% number of subsystems
n = length(A);

% to speed up the algorithm a bit
cache{1} = 1;
for k=1:n
  dim(k) = length(A{k});
  D(k) = size(A{k}{1}, 1);

  % initialize matrix product cache
  cache{k+1} = cache{k} * A{k}{1};
end


% contractions
ind = zeros(1, n); % subsystem indices
for k=1:prod(dim)
  x(k) = trace(cache{n+1});
  
  % increment subsys indices
  p = n;
  while (p >= 1)
    ind(p) = ind(p)+1;
    if (ind(p) >= dim(p))
      ind(p) = 0;
      p = p-1; 
    else
      % recompute the changed subproducts
      for j=p:n
        cache{j+1} = cache{j} * A{j}{ind(j)+1};
      end
      break
    end
  end
end

s = normalize(state(x, dim));
