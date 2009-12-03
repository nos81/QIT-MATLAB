function s = tensor(varargin)
% STATE/TENSOR  Tensor product of states.
%  s = tensor(s1, s2, ...)
%
%  Returns the tensor product state of states s1, s2, ...

% Ville Bergholm 2009


s = varargin{1};
s.dim = [];
pure = true;

% concatenate dim vectors
for k = 1:nargin
  s.dim = cat(2, s.dim, varargin{k}.dim);

  % if all states are pure, keep the result state pure
  if (size(varargin{k}.data, 2) ~= 1)
    pure = false;
  end
end


s.data = 1;
for k = 1:nargin
  if (~pure && (size(varargin{k}.data, 2) == 1))
    % vector into matrix
    s.data = kron(s.data, varargin{k}.data*varargin{k}.data');
  else
    s.data = kron(s.data, varargin{k}.data);
  end
end
