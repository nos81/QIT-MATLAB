function s = tensor(varargin)
% TENSOR  Tensor product of states.
%  s = tensor(s1, s2, ...)
%
%  Returns the tensor product state of states s1, s2, ...

% Ville Bergholm 2009-2010


% if all states are pure, keep the result state pure
pure = true;
for k = 1:nargin
  if (size(varargin{k}.data, 2) ~= 1)
    pure = false;
    break;
  end
end

% otherwise convert all states to state ops before tensoring
if (~pure)
  for k = 1:nargin
    varargin{k} = to_op(varargin{k});
  end
end

s = tensor@lmap(varargin{:});

s = remove_singletons(s);
