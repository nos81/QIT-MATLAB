function I = LU(rho, k, perms)
% Local unitary polynomial invariants of quantum states.
%  I = LU(rho, k, perms)
%
%  Computes the permutation invariant I_{k; p1, p2, ..., pn} for the state rho.
%  perms is a cell vector containing n permutation vectors, one for
%  each subsystem in the state.
%  NOTE: Full permutation vectors, NOT cycles!
%
%  k is the order/degree of the invariant (number of copies of rho in the
%  corresponding diagram). Each pj must thus be a k-permutation.
%
%  Example: I_{3; (123),(12)}(rho) = LU_inv(rho, 3, {[2 3 1], [2 1 3]})
%
%  This function can be very inefficient for some invariants, since
%  it does no partial traces etc. which might simplify the calculation.

% shortcut: permutation [a b] means swap a with b...

% Ville Bergholm 2011-2014


n = length(perms);
if n ~= rho.subsystems()
    error('Need one permutation per subsystem.')
end

for j = 1:n
  perms{j} = perm_convert(perms{j}, k);
end

% splice k sequential copies of the entire system into k copies of each subsystem
s = reshape(1:n*k, [k, n]).';
s = s(:).';

% permute the k copies of each subsystem
temp = kron(k*(0:n-1), ones(1, k));
p = cell2mat(perms) + temp;

% Permutations: a*b = a(b), x = y * z^{-1}  <=>  x * z = x(z) = y.
s_inv(s) = 1:n*k;
total = s_inv(p(s)); % total = s^{-1} * p * s

% TODO this could be done much more efficiently
temp = lmap(tensorpow(rho.to_op(), k));
I = trace(reorder(temp, {total, []}));
end


function ret = perm_convert(p, k)
% interpret shortcut permutations

id = 1:k; % identity permutation

% troublesome case (not a swap!)
if k == 2 && isequal(p, [1 2])
  p = [];
end

switch length(p)
  case 0
    % [], identity
    ret = id;

  case 2
    % swap two subsystems
    ret = id;
    ret(p(1)) = p(2);
    ret(p(2)) = p(1);
  
  otherwise
    % full permutation
    if length(setxor(id, p)) ~= 0
        error('Invalid permutation.');
    end
    ret = p;
end
end
