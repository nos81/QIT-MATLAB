function I = LU_inv(rho, k, perms)
% Local unitary polynomial invariants of quantum states.
%  I = LU_inv(rho, k, perms)
%
%  Computes the permutation invariant I_{k; p1, p2, ..., pn} for the state rho.
%  perms is a cell vector containing n permutation vectors.
%
%  Example: I_{3; (123),(12)}(rho) = LU_inv(rho, 3, {[2 3 1], [2 1 3]})
%
%  This function can be very inefficient for some invariants, since
%  it does no partial traces etc. which might simplify the calculation.
    
% Ville Bergholm 2011


n = length(perms);
d = dims(rho);
if n ~= length(d)
    error('Need one permutation per subsystem.')
end

% convert [] to identity permutation
for j = 1:n
    if length(perms{j}) == 0
        perms{j} = 1:k;
    end
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
I = trace(reorder(tensor_pow(rho, k), {total, []}));
end



function ret = tensor_pow(rho, n)
% Returns \rho^{\otimes n}.

rho = to_op(rho);
ret = lmap(rho);
for k=2:n
    ret = tensor(ret, rho);
end
end
