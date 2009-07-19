classdef state
% STATE  Class for quantum states.
%
%  Describes a discrete composite quantum system consisting of subsystems
%  whose dimensions are given in the vector dim (big-endian ordering).
%  Can handle both pure and mixed states.

% Ville Bergholm 2008-2009

  properties
    data % state vector/operator
    dim  % subsystem dimensions vector
  end

  methods
    function out = state(s, dim)
    % constructor
    %  x = state(s [, dim]);
    %  x = state(y);                % copy constructor
    %  x = state(y, [2 2]);         % copy constructor, reinterpret system as two qubits
    %  x = state("00101");          % standard basis ket |00101> in a five-qubit system
    %  x = state("02", [2 3]);      % standard basis ket |02> in a qubit+qutrit system
    %  x = state(k, [2 3]);         % standard basis ket |k> in a qubit+qutrit system, k must be an integer scalar
    %  x = state(rand(4,1));        % ket, dim = 4
    %  x = state(rand(4,1), [2 2]); % ket, two qubits
    %  x = state(rand(4));          % state operator, dim = 4
    %  x = state(rand(6), [3 2]);   % state operator, qutrit+qubit
    %  x = state("GHZ", [2 2 2]);   % three-qubit GHZ state
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
        ind = s(1) - '0';
        for k=2:n
          ind = (ind * dim(k)) + (s(k) - '0');
        end
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
        s = zeros(prod(dim), 1);
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

    % error checking
    if (size(s, 1) ~= prod(dim))
      error('Dimension of the system does not match the combined dimension of the subsystems.')
    end

    out.data = s;
    out.dim = dim; % big-endian ordering
    end

    function x = subsref(s, index)
    % STATE/SUBSREF  Direct access to the data members.
      switch (index.type)
	case '.'
	  switch (index.subs)
	    case 'data'
	      x = s.data;
	    case 'dim'
	      x = s.dim;
  	    otherwise
	      error('Unknown state property.')
	  end

	otherwise
	  error('State class cannot be indexed with that operator.')
      end
    end

    function sys = clean_selection(s, sys)
    % STATE/CLEAN_SELECTION  Internal helper, makes a subsystem set unique and sorted.
      sys = intersect(1:length(s.dim), sys);
    end
  end
end
