classdef state < lmap
% STATE  Class for quantum states.
%
%  Describes a discrete composite quantum system consisting of subsystems
%  whose dimensions are given in the vector dim (big-endian ordering).
%  Can handle both pure and mixed states.

%  State class instances are special cases of lmaps. They have exactly two indices.
%  TODO If the dimension of the first  index is 1, it is a bra.
%  If the dimension of the second index is 1, it is a ket.
%  If neither of the above is true, both indices must have equal dimensions and the
%  object represents a state operator.

% Ville Bergholm 2008-2010

  methods
    function out = state(s, dim)
    % constructor
    %  x = state(s [, dim]);
    %  x = state(y);                % copy constructor
    %  x = state(y, [2 2]);         % copy constructor, reinterpret system as two qubits
    %  x = state('00101');          % standard basis ket |00101> in a five-qubit system
    %  x = state('02', [2 3]);      % standard basis ket |02> in a qubit+qutrit system
    %  x = state(k, [2 3]);         % standard basis ket |k> in a qubit+qutrit system, k must be an integer scalar
    %  x = state(rand(4,1));        % ket, dim = 4
    %  x = state(rand(4,1), [2 2]); % ket, two qubits
    %  x = state(rand(4));          % state operator, dim = 4
    %  x = state(rand_positive(6), [3 2]); % state operator, qutrit+qubit
    %  x = state('GHZ', [2 2 2]);   % three-qubit GHZ state
    %
    %  The currently supported named states are GHZ (Greenberger-Horne-Zeilinger) and W.

    if (nargin == 0)
      error('No arguments given.');

    elseif (nargin > 2)
      error('Too many arguments.');

    elseif (isa(s, 'state'))
      % copy constructor
      if (nargin == 1)
        dim = s.dim;  % copy also dimensions
      end
      s = s.data;

    elseif (isa(s, 'lmap'))
      % state from lmap
      if (nargin == 1)
        dim = s.dim;  % copy also dimensions
	    dim = dim{1};
      end
      s = s.data;

    elseif (ischar(s))
      % string
      if (isletter(s(1)))
        % named state
        if (nargin == 1)
          dim = [2 2 2];
        end

        name = lower(s);
        n = length(dim);
        s = zeros(prod(dim), 1);
        switch name
          case {'bell1', 'bell2', 'bell3', 'bell4'}
            Q_Bell = [1 0 0 i; 0 i 1 0; 0 i -1 0; 1 0 0 -i] / sqrt(2);
            dim = [2 2];
            s = Q_Bell(:, name(5)-'0');
            
          case 'ghz'
            s(1) = 1; s(end) = 1;

          case 'w'
            ind = 1;
            for k=n:-1:1
              s(ind*(dim(k) - 1) + 1) = 1; % MATLAB indexing starts at 1
              ind = ind*dim(k);
            end

          otherwise
            error('Unknown named state ''%s''.', name)
        end
        s = s/norm(s); % normalize
      
      else
        % number string defining a standard basis ket
        if (nargin == 1)
          n = length(s); % number of units
          dim = 2*ones(1, n); % assume units are qubits
        end

        % calculate the linear index
        n = length(dim);
        s = s - '0';
        if (any(s >= dim))
          error('Invalid basis ket.')
        end
        muls = fliplr(circshift(cumprod(fliplr(dim)), [0 1]));
        muls(end) = 1;
        ind = muls*s.';
        s = zeros(prod(dim), 1);
        s(ind+1) = 1; % MATLAB indexing starts at 1
      end

    elseif (isnumeric(s))
      if (isscalar(s))
        % integer represented using the computational standard basis
        if (nargin == 1)
          error('Need system dimension.')
        end

        ind = s;
        temp = prod(dim); % total number of states
        if (ind >= temp)
          error('Invalid basis ket.')
        end
        s = zeros(temp, 1);
        s(ind+1) = 1; % MATLAB indexing starts at 1

      else
        % state vector or matrix
        sss = size(s);
        if (length(sss) > 2)
          error('State must be given as a state vector or a state operator.')
        end
        if (sss(1) == 1)
          s = s.'; % row vector into column vector
        elseif (sss(2) > 1 && sss(1) ~= sss(2))
          error('State operator matrix must be square.')
        end
      end

      if (nargin == 1)
        dim = size(s, 1); % no dim given, infer from s
      end
    end

    % state vector or operator?
    if (size(s, 2) ~= 1)
      dim = {dim, dim};
    else
      dim = {dim, 1};
    end

    % call the lmap constructor
    out = out@lmap(s, dim);
    end

    function x = subsref(s, index)
    % SUBSREF  Direct access to the data members.
      switch (index.type)
        case '.'
          switch (index.subs)
            case 'data'
              x = s.data;
            case 'dim'
              x = s.dim{1}; % since for a state the indices must be identical or singletons
            otherwise
              error('Unknown state property.')
          end

        otherwise
          error('State class cannot be indexed with that operator.')
      end
    end

    function d = dims(s)
      d = s.dim{1}; % dims of the other index must be equal (or 1)
    end

    function ret = equal_dims(s, t)
      ret = isequal(dims(s), dims(t));
    end

    function ret = is_ket(s)
      ret = (size(s.data, 2) == 1);
    end
  end
end
