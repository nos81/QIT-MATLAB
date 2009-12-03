function a = bloch_vector(s)
% STATE/BLOCH_VECTOR  Generalized Bloch vector corresponding to the state.
%  a = bloch_vector(s)
%
%  Returns the Bloch vector corresponding to state s.
%  The vector is defined in terms of the tensor basis
%  corresponding to s.dim.
%
%  For proper states norm(a) <= sqrt(prod(s.dim)-1) (e.g. for qubits norm(a) <= 1).

% Ville Bergholm 2009


d = prod(s.dim);
n = d^2-1;
G = tensorbasis(s.dim);

for k=1:n
  a(k) = ev(s, G{k});
end
a = a'*sqrt(d); % to match the usual Bloch vector normalization
