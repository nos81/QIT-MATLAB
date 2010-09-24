classdef lmap
% LMAP  Class for linear maps between lists of finite-dimensional Hilbert spaces.
%
%  Contains both the tensor and the dimensional information.

% Ville Bergholm 2008-2010

  properties
    data % array of tensor elements
    dim  % cell vector of subsystem dimension row vectors, one for each index
         % the dimensions are in big-endian order, the indices in standard matrix order
         % right now we only support one or two indices, i.e. vectors and matrices
  end

  methods
    function out = lmap(s, dim)
    % constructor
    %  x = lmap(s [, dim]);
    %    = lmap(y);                % (y is an lmap) copy constructor
    %    = lmap(y, dim);           % (y is an lmap) copy constructor, reinterpret dimensions
    %    = lmap(k, dim);           % standard basis ket |k> or bra <k|. k must be an integer scalar
    %    = lmap(rand(4,1));        % ket, dim = {4, 1}
    %    = lmap(rand(4,1), [2 2]); % ket, dim = {[2 2], 1}
    %    = lmap(rand(4));          % operator, dim = {4, 4}
    %    = lmap(rand(6), [3 2]);   % operator

    if (nargin == 0)
      error('No arguments given.');
    elseif (nargin > 2)
      error('Too many arguments.');
    end

    if (nargin == 2)
      % cover for lazy users
      if (isnumeric(dim))
        dim = {dim};
      end
      if (length(dim) < 2)
        dim{2} = 1;
      end
    end

    if (isa(s, 'lmap'))
      % copy constructor
      if (nargin == 1)
        dim = s.dim;  % copy also dimensions
      end
      s = s.data;

    elseif (isnumeric(s))
      if (isscalar(s))
        % integer denoting a computational basis ket
        if (nargin == 1)
          error('Need system dimensions.')
        end

        % NOTE only works because we have max. two indices
        if (isequal(dim{1}, 1))
          q = 2;
        elseif (isequal(dim{2}, 1))
          q = 1;
        else
          error('Construction by ket/bra number only for vectors.')
        end

        ind = s;
        d = [prod(dim{1}), prod(dim{2})]; % total number of states

        if (ind >= d(q))
          error('Invalid basis ket/bra.')
        end
        s = zeros(d);
        s(ind+1) = 1; % MATLAB indexing starts at 1

      else
        % full vector or matrix
        sss = size(s);
        if (length(sss) > 2)
          error('For now only vectors and matrices.')
        end

        if (nargin == 1)
          dim = num2cell(size(s)); % no dim given, infer from s
        end
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
      error('not enough dims given');
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
      sys = intersect(1:length(s.dim), sys);
    end
  end
end
