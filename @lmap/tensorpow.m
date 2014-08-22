function ret = tensorpow(U, n)
% TENSORPOW  Tensor power of the lmap.
%  x = tensorpow(U, n)
%
%  Returns the tensor power U^{\otimes n}.

% Ville Bergholm 2014


if n < 1
    error('Only positive integer tensor powers are allowed.')
end

dd = order(U);
ret = U;
for k = 2:n
    for ind = 1:dd
        ret.dim{ind} = cat(2, ret.dim{ind}, U.dim{ind}); % concatenate dimension lists
    end

    % kronecker product of the data
    ret.data = kron(ret.data, U.data);
end

ret = ret.remove_singletons();
