classdef lmap
% LMAP  Class for multilinear maps between tensor products of finite-dimensional Hilbert spaces.
%
%  Contains both the tensor and the dimensional information.

%  TODO Another possible interpretation of lmap would be to only consider order-2 lmaps and
%  treat each subsystem as an index, with the subsystems within dim{1} and dim{2}
%  corresponding to contravariant and covariant indices, respectively?

% Ville Bergholm 2008-2011

  properties
    data % array of tensor elements
    dim  % cell vector of subsystem dimension row vectors, one for each index
         % the dimensions are in big-endian order, the indices in standard matrix order
  end

  methods
    function out = lmap(s, dim)
    % constructor
    %  x = lmap(s [, dim]);             % dim == {out_dims, in_dims}
    %    = lmap(y);                     % (y is an lmap) copy constructor
    %    = lmap(y, dim);                % (y is an lmap) copy constructor, reinterpret dimensions
    %    = lmap(rand(4,1));             % ket, 1 -> 4
    %    = lmap(rand(4,1), {[2 2], 1}); % ket, 1 -> [2 2]
    %    = lmap(rand(4));               % operator, 4 -> 4
    %    = lmap(rand(6), {6, [3 2]});   % operator, [3 2] -> 6
 
    if nargin == 0
      error('No arguments given.');
    elseif nargin > 2
      error('Too many arguments.');
    end

    if nargin == 2
      % cover for lazy users, convert a dim vector into a cell vector with one entry
      if isnumeric(dim)
        dim = {dim};
      end
    end

    if isa(s, 'lmap')
      % copy constructor
      if nargin == 1
        dim = s.dim; % copy also dimensions
      end
      s = s.data;

    elseif isnumeric(s)
      % full tensor
      if nargin == 1
        dim = num2cell(size(s)); % no dim given, infer from s
      end
    end

    % error checking
    n = length(dim);
    for k=1:n
      if size(s, k) ~= prod(dim{k})
        error('Dimensions of index %d do not match the combined dimensions of the subsystems.', k)
      end
    end
    if size(s, k+1) ~= 1
      error('The dimension cell vector does not have an entry for each tensor index.');
    end

    out.data = s;
    out.dim = dim; % big-endian ordering
    end


    function s = remove_singletons(s)
    % REMOVE_SINGLETONS  Eliminate unnecessary singleton dimensions.

      for k = 1:length(s.dim)
	s.dim{k}(find(s.dim{k} == 1)) = []; % remove all ones
	if isempty(s.dim{k})
	  s.dim{k} = 1; % restore a single one if necessary
	end
      end
    end


    function n = order(s)
    % ORDER  Tensor order of the lmap.
      n = length(s.dim);
    end


    function ret = is_compatible(s, t)
    % IS_COMPATIBLE  True iff s and t have same order and equal dimensions and can thus be added.

      n = order(s);
      m = order(t);

      if n ~= m
        error('The orders of the lmaps do not match.')
      end

      for k = 1:n
        if ~isequal(s.dim{k}, t.dim{k})
          error('The dimensions of the index %d of the lmaps do not match.', k)
        end
      end
      ret = true;
    end


    function ret = is_concatenable(s, t)
    % IS_CONCATENABLE  True iff s * t is a meaningful expression.
    %  ret = is_concatenable(s, t)
    %  
    %  s and t can be concatenated (multiplied) if they are both
    %  order-2 lmaps and the input dimensions of s match the output
    %  dimensions of t.

      if order(s) ~= 2 || order(t) ~= 2
        error('Concatenation only defined for order-2 lmaps.')
      end

      ret = isequal(s.dim{2}, t.dim{1});
      if ~ret
          error('The input and output dimensions do not match.')
      end
    end
    
    
    function ret = is_ket(s)
    % IS_KET  Returns true iff the lmap is a ket vector.
        ret = (size(s.data, 2) == 1);
    end
    
    
    function x = norm(s)
    % NORM  Matrix norm of the lmap.
    %  x = norm(s)
    %
    %  Returns the 2-norm of the lmap s.
      x = norm(s.data);
    end


    function s = conj(s)
    % CONJ  Complex conjugate of the lmap.
    %  q = conj(s)
    %
    %  Returns the lmap complex conjugated in the computational basis.

        s.data = conj(s.data);
    end

    function s = real(s)
    % REAL  Real part of the lmap.
    %  q = real(s)
    %
    %  Returns the real part of the lmap in the computational basis.

        s.data = real(s.data);
    end

    function s = imag(s)
    % IMAG  Imaginary part of the lmap.
    %  q = imag(s)
    %
    %  Returns the imaginary part of the lmap in the computational basis.

        s.data = imag(s.data);
    end


    function x = trace(s)
    % TRACE  Trace of the lmap.
    %  x = trace(s)
    %
    %  Returns the trace of the lmap s.
    %  The trace is only properly defined if s.dim{1} == s.dim{2}.

      if ~isequal(s.dim{1}, s.dim{2})
        error('Trace only defined for endomorphisms ("square matrices").')
      end
    
      x = trace(s.data);
    end
  end
end
