function s = tensor(varargin)
% TENSOR  Tensor product of lmaps.
%  s = tensor(s1, s2, ...)
%
%  Returns the tensor product of lmaps s1, s2, ...

% Ville Bergholm 2009-2010


s = varargin{1};

for k = 2:nargin
  n = order(varargin{k});

  if (order(s) < n)
    s.dim{n} = []; % all missing indices filled in with []
  end

  for ind = 1:n
    s.dim{ind} = cat(2, s.dim{ind}, varargin{k}.dim{ind}); % concatenate dimension lists
  end

  % kronecker product of the data
  s.data = kron(s.data, varargin{k}.data);
end

%s = remove_singletons(s);
