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
    %  x = lmap(s [, dim]);
    %    = lmap(y);                % (y is an lmap) copy constructor
    %    = lmap(y, dim);           % (y is an lmap) copy constructor, reinterpret dimensions
    %    = lmap(k, dim);           % standard basis ket |k> or bra <k|. k must be an integer scalar
    %    = lmap(rand(4,1));        % ket, dim = {4, 1}
    %    = lmap(rand(4,1), [2 2]); % ket, dim = {[2 2]}
    %    = lmap(rand(4));          % operator, dim = {4, 4}
    %    = lmap(rand(6), [3 2]);   % operator

    if (nargin == 0)
      error('No arguments given.');
    elseif (nargin > 2)
      error('Too many arguments.');
    end

    if (nargin == 2)
      % cover for lazy users, convert a dim vector into a cell vector with one entry
      if (isnumeric(dim))
        dim = {dim};
      end
    end

    if (isa(s, 'lmap'))
      % copy constructor
      if (nargin == 1)
        dim = s.dim; % copy also dimensions
      end
      s = s.data;

    elseif (isnumeric(s))
      % full tensor
      if (nargin == 1)
        dim = num2cell(size(s)); % no dim given, infer from s
      end
    end

    % error checking
    n = length(dim);
    for k=1:n
      if (size(s, k) ~= prod(dim{k}))
        error('Dimensions of index %d do not match the combined dimensions of the subsystems.', k)
      end
    end
    if (size(s, k+1) ~= 1)
      error('The dimension cell vector does not have an entry for each tensor index.');
    end

    out.data = s;
    out.dim = dim; % big-endian ordering
    end


    function x = subsref(s, index)
    % SUBSREF  Direct access to the data members.
      switch (index.type)
	case '.'
	  switch (index.subs)
	    case 'data'
	      x = s.data;
	    case 'dim'
	      x = s.dim;
  	    otherwise
	      error('Unknown lmap property.')
	  end

	otherwise
	  error('Lmap class cannot be indexed with that operator.')
      end
    end


    function sys = clean_selection(s, sys)
    % CLEAN_SELECTION  Internal helper, makes a subsystem set unique and sorted.
      sys = intersect(1:length(s.dim{1}), sys);
    end


    function s = remove_singletons(s)
    % REMOVE_SINGLETONS  Eliminate unnecessary singleton dimensions.

      for k = 1:length(s.dim)
	s.dim{k}(find(s.dim{k} == 1)) = []; % remove all ones
	if (isempty(s.dim{k}))
	  s.dim{k} = 1; % restore a single one if necessary
	end
      end
    end


    function n = order(s)
    % ORDER  Tensor order of the lmap.
      n = length(s.dim);
    end


    function ret = is_compatible(s, t)
    % IS_COMPATIBLE  True iff s and t have same order and equal dimensions.

      n = order(s);
      m = order(t);

      if (n ~= m)
        error('The orders of the lmaps do not match.')
      end

      for k = 1:n
        if ~isequal(s.dim{k}, t.dim{k})
          error('The dimensions of the index %d of the lmaps do not match.', k)
        end
      end
      ret = true;
    end

    
    function x = norm(s)
    % NORM  Norm of the lmap.
    %  x = norm(s)
    %
    %  Returns the 2-norm of the lmap s.

      x = norm(s.data);
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
