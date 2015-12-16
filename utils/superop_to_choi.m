function [C] = superop_to_choi(L)
% SUPEROP_TO_CHOI  Convert a Liouvillian superoperator to a Choi matrix.
%  [C] = superop_to_choi(L)
%
% Given a Liouvillian superoperator L operating on vectorized state
% operators, L: A \otimes A \to B \otimes B, returns the
% corresponding Choi matrix C: A \otimes B \to A \otimes B.

% Ville Bergholm 2015


% input and output Hilbert space dimensions
s = size(L);
dim = sqrt(s);
A = dim(2);
B = dim(1);

% reshape uses little-endian ordering, whereas we use big-endian
temp = reshape(L, [B, B, A, A]);
% swap the more significant output and less significant input indices
temp = permute(temp, [1, 3, 2, 4]);
% back to a flat matrix
C = reshape(temp, [A*B, A*B]);
