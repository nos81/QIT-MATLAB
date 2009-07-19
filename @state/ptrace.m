function s = ptrace(s, sys)
% STATE/PTRACE  Partial trace.
%  x = ptrace(s, sys)
%  x = ptrace(s, [2 3]) % trace over subsystems 2 and 3 of state s
%
%  Returns the partial trace of the state s
%  wrt. the subsystems listed in the vector sys.

% Ville Bergholm 2008-2009


% number of systems
n = length(s.dim);

% we trace over the subsystems in order, starting from the first one
sys = clean_selection(s, sys);

% big-endian ordering is more natural for users, but little-endian more convenient for calculations
d = fliplr(s.dim); % into little-endian: now dim(1) is the dimension of the last system (least significant digit)
flipped_sys = n+1-sys; % these are indices to dim vector, flip them too, now d(flipped_sys) makes sense

if (size(s.data, 2) == 1)
  s.data = s.data*s.data'; % vector into matrix
end

% partial trace over single system s, performed for every s in sys
for j=flipped_sys
  muls = circshift(cumprod(d), [0 1]); % index multipliers
  muls(end+1) = muls(1); % convenience
  muls(1) = 1; % now muls == [1, d_n, d_n*d_{n-1}, ... ]

  % build the index "stencil"
  inds = 1; % matlab indexing starts at 1
  for k=1:n
    if (k == j)
      continue
    end
    inds = tensorsum([0:muls(k):(muls(k+1)-1)], inds);
  end

  stride = muls(j); % stride for moving the stencil while summing
  res = zeros(length(inds)); % result

  for k=0:(d(j)-1)
    temp = inds + stride*k;
    res = res + s.data(temp, temp);
  end

  s.data = res; % replace data
  d(j) = 1;  % remove traced-over dimension.
end

s.dim = s.dim(setdiff(1:n, sys)); % remove traced-over dimensions for good

return



function c = tensorsum(a,b)

%c = log(kron(exp(a), exp(b))); % a perverted way of doing it, the exp overflows...

c = [];
for k=1:length(a)
  c = [c, a(k)+b];
end

return
