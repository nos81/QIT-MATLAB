function s = coherent_state(alpha, n)
% COHERENT_STATE  Coherent states of a harmonic oscillator.
%  s = coherent_state(alpha, n)
%
%  Returns the n-dimensional approximation to the
%  coherent state \ket{\alpha}.

% Ville Bergholm 2010


ket = zeros(n, 1);
for k=0:n-1
  ket(k+1) = alpha^k/sqrt(factorial(k));
end
%ket = ket*exp(-abs(alpha)^2/2);
s = normalize(state(ket, n));

return

s = exp(-abs(alpha)^2/2) * u_propagate(state(0, n), expm(alpha*ho.ladder(n)'));
s = u_propagate(state(0, n), ho.displace(alpha, n));

% normalization!?
