function s = tensor(varargin)
% TENSOR  Tensor product of lmaps.
%  s = tensor(s1, s2, ...)
%
%  Returns the tensor product of lmaps s1, s2, ...

% Ville Bergholm 2009-2010


s = varargin{1};

% concatenate dim vectors
for k = 2:nargin
  n = length(varargin{k}.dim);

  if (length(s.dim) < n)
    s.dim{n} = []; % all missing indices filled in with []
  end

  for ind = 1:n
    %if (isequal(s.dim{ind}, 1))
    %  s.dim{ind} = []; % get rid of singletons
    %end
    s.dim{ind} = cat(2, s.dim{ind}, varargin{k}.dim{ind}); % concatenate dimension lists
  end

  s.data = kron(s.data, varargin{k}.data);
end

% eliminate extra singleton dimensions
for k=1:length(s.dim)
  s.dim{k}(find(s.dim{k} == 1)) = [];
  if (isempty(s.dim{k}))
    s.dim{k} = 1; % restore one
  end
end
