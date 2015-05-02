function C = copydot(n_in, n_out, d)
% COPYDOT  Copy dot operators.
%  C = copydot(n_in, n_out, d)
%
%  Returns the matrix representation for the n_in -> n_out copy dot
%  with subsystem dimension d.
    
%! See reference BB2011.
% Ville Bergholm 2010-2014


if nargin < 3
    d = 2;  % qubits by default
end

C = sparse(d^n_out, d^n_in);

% compute the strides by summing up 1+d+d^2+...+d^(n-1)
stride_in  = (d^n_in -1) / (d-1);
stride_out = (d^n_out-1) / (d-1);

% loop over the sum
for k=0:d-1
    C(1 +k*stride_out, 1 +k*stride_in) = 1;
end

C = lmap(C, {d * ones(1, n_out), d * ones(1, n_in)});
