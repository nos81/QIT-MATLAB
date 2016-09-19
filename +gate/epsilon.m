function U = epsilon(n)
% EPSILON  Epsilon tensor.
%  U = epsilon(d)
%
%  Returns the fully antisymmetric Levi-Civita symbol in n dimensions,
%  a tensor with n n-dimensional subsystems.

% Ville Bergholm 2016

D = n^n;
dim = n * ones(1,n);
U = sparse(D,1);

p = perms(1:n);
% loop through all permutations of 1:n
for k=1:size(p,1)
    temp = p(k,:);
    ind = ravel_index(temp, dim);
    U(ind) = sign(temp);
end
U = lmap(U, {dim, 1});
end

function s = sign(p)
% sign of the permutation p by counting inversions
    n = length(p);
    s = 0;
    for k=1:n
        for j=k+1:n
            s = s +(p(k)>p(j));
        end
    end
    s = mod(s, 2);
    s = (-1)^s;
end
