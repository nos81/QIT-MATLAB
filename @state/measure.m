function [p, res, s] = measure(s, M, discard)
% MEASURE  Quantum measurement.
%  [p, res, s]
%    = measure(s)                % measure the entire system projectively
%    = measure(s, [1 4])         % measure subsystems 1 and 4 projectively
%    = measure(s, {M1, M2, ...}) % perform a general measurement
%    = measure(s, A)             % measure a Hermitian observable A
%
%  Performs a quantum measurement on the state s.
%
%  If only s is given, a full projective measurement in the
%  computational basis is performed.
%
%  If a vector of subsystems is given as the second parameter, only
%  those subsystems are measured, projectively, in the
%  computational basis.
%
%  A general measurement may be performed by giving a complete set
%  of measurement operators {M1, M2, ...} as the second parameter.
%
%  Finally, if the second parameter is a Hermitian matrix A, the
%  corresponding observable is measured. In this case the second
%  column of p contains the eigenvalue of A corresponding to each
%  measurement result.
%
%  p = measure(...) returns the vector p, where p(k) is the probability of
%  obtaining result k in the measurement. For a projective measurement
%  in the computational basis this corresponds to the ket |k-1>.
%
%  [p, res] = measure(...) additionally returns the index of the result of the
%  measurement, res, chosen at random following the probability distribution p.
% 
%  [p, res, s] = measure(...) additionally gives s, the collapsed state
%  corresponding to the measurement result res.

% Ville Bergholm 2009-2010


d = s.dim{1};

if (nargin <= 1)
  % full measurement in the computational basis
  p = prob(s); % probabilities 

  if (nargout >= 2)
    res = rand_measure(p);
    if (nargout >= 3)
      s = state(res-1, d); % collapsed state
    end
  end
  return;

elseif (isnumeric(M))
 if (isvector(M))
  % measure a set of subsystems in the computational basis
  sys = clean_selection(s, M);

  % dimensions of selected subsystems and identity ops between them
  % TODO sequential measured subsystems could be concatenated as well
  q = length(sys);
  ppp = 1; % first sys not yet included
  for k=1:q
    dims(2*k-1) = prod(d(ppp:sys(k)-1)); % identity
    dims(2*k) = d(sys(k)); % selected subsys
    ppp = sys(k)+1;
  end
  dims(2*q+1) = prod(d(ppp:end)); % last identity

  % big-endian ordering is more natural for users, but little-endian more convenient for calculations
  muls = fliplr(circshift(cumprod(fliplr(d(sys))), [0 1]));
  m = muls(end); % number of possible results
  muls(end) = 1; % now muls == [..., d_s{q-1}*d_s{q}, d_s{q}, 1]

  % sum the probabilities
  born = prob(s);
  for j=1:m
    stencil = build_stencil(j, q, dims, muls);
    p(j) = stencil*born; % inner product
  end

  if (nargout >= 2)
    res = rand_measure(p);
    if (nargout >= 3)
      R = build_stencil(res, q, dims, muls); % each projector is diagonal, hence we only store the diagonal

      if (nargin == 3 && discard)
        % discard the measured subsystems from s
        
        d(sys) = [];
        keep = find(R);  % indices of elements to keep
        
        if is_ket(s)
          % state vector
          s.data = s.data(keep) / sqrt(p(res)); % collapsed state
        else
          % state operator
          s.data = s.data(keep, keep) / p(res); % collapsed state
        end
        s = state(s, d);
      else
        if is_ket(s)
          % state vector
          s.data = R.' .* s.data / sqrt(p(res)); % collapsed state
        else
          % state operator
          s.data = (R.'*R) .* s.data / p(res); % collapsed state, HACK
        end
      end
    end
  end

 else % M is a matrix
  % measure the given Hermitian observable

  [a, P] = spectral_decomposition(M);
     
  m = length(a); % number of possible results
  p = zeros(m, 2);
  for k=1:m
    p(k, 1) = ev(s, P{k});
  end
  p(:,2) = a; % also return the corresponding results

  if (nargout >= 2)
    res = rand_measure(p);
    if (nargout >= 3)
      ppp = P{res};

      if is_ket(s)
        % state vector
        s.data = ppp * s.data / sqrt(p(res)); % collapsed state
      else
        % state operator
        s.data = ppp * s.data * ppp / p(res); % collapsed state
      end
    end
  end
 end
elseif (iscell(M))
  % otherwise use set M of measurement operators (assumed complete!)
  m = length(M);

  % probabilities
  if is_ket(s)
    % state vector
    for k=1:m
      p(k) = s.data' * M{k}' * M{k} * s.data;
    end
    
    if (nargout >= 2)
      res = rand_measure(p);
      if (nargout >= 3)
        s.data = M{res} * s.data / sqrt(p(res)); % collapsed state
      end
    end

  else
    % state operator
    for k=1:m
      p(k) = trace(M{k}' * M{k} * s.data); % TODO wasteful
    end
    
    if (nargout >= 2)
      res = rand_measure(p);
      if (nargout >= 3)
        s.data = M{res} * s.data * M{res}' / p(res); % collapsed state
      end
    end
  end

else
  error('Unknown input type.')
end


end


function stencil = build_stencil(j, q, dims, muls)
% build projector to state j (diagonal because we project into the computational basis)

stencil = ones(1, dims(1)); % first identity
for k=1:q
  temp = sparse(1, dims(2*k));
  temp(mod(floor((j-1)/muls(k)), dims(2*k))+1) = 1; % projector
  stencil = kron(kron(stencil, temp), ones(1, dims(2*k+1))); % identity
end
end

function r = rand_measure(p)
% random measurement

  r = find(cumsum(p) >= rand, 1, 'first');
end
