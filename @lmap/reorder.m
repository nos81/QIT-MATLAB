function [s] = reorder(s, perm)
% REORDER  Change the relative order of subsystems in an lmap.
%  [S] = reorder(T, perm);
%    reorder(T, {[], [3 2 1]}); % ignore first index, reverse the order of subsystems in the second
%    reorder(T, {[2 5]});       % swap subsystems 2 and 5 in the first index
%
%  Reorders the subsystems of the lmap T according to permutation
%  cell vector perm == {p_1, p_2, ...}.
%  perm contains one permutation vector p_i for each index i of T.
%  If there are more indices than permutation vectors, the extra indices remain untouched.
%
%  The contents of each permutation vector p_i may be either empty (nothing is done to that index), 
%  exactly two subsystem numbers (to be swapped), or a full permutation of subsystem numbers.

% Ville Bergholm 2009-2010


% number of indices
n_ind = order(s);

if (length(perm) < n_ind)
  perm{n_ind} = []; % all missing permutations filled in with []
end

total_d = [];
total_perm = [];
last_used_index = 0;

% loop over indices
for k = 1:n_ind

  this_perm = perm{k};  % requested subsystem permutation for this index
  this_dim  = s.dim{k}; % subsystem dims for this index

  % number of subsystems
  n = length(this_dim);

  % total dimension
  dd(k) = prod(this_dim);

  temp = 1:n;
  switch (length(this_perm))
    case 0
      % [], no change
      % let the dimensions vector be, lump all subsystems in this index into one
      this_dim = dd(k);
      this_perm = 1;
      n = 1;

    case 2
      % swap two subsystems
      temp(this_perm(1)) = this_perm(2);
      temp(this_perm(2)) = this_perm(1);
      this_perm = temp;

      s.dim{k} = this_dim(this_perm); % reorder the dimensions vector

    otherwise
      % full permutation
      if (length(setxor(temp, this_perm)) ~= 0)
        error('Invalid permutation for index %d.', k);
      end

      s.dim{k} = this_dim(this_perm); % reorder the dimensions vector
  end

  % big-endian ordering is more natural for users, but Matlab funcs
  % prefer little-endian, so we reverse it
  total_d    = [total_d,    fliplr(this_dim)];
  total_perm = [total_perm, last_used_index + fliplr(n+1 -this_perm)];
  last_used_index = last_used_index + n;
end


T = full(s.data); % FIXME Matlab 7.6 reshape is broken (sparse matrices cause problems with singleton dimensions)

% tensor into another tensor which has one index per subsystem, permute dimensions, back into a tensor with the original number of indices
s.data = reshape(permute(reshape(T, total_d), total_perm), [dd 1]);
