function s = coherent_state(alpha, n)
% COHERENT_STATE  Coherent states of a harmonic oscillator.
%  s = coherent_state(alpha, n)
%
%  Returns the n-dimensional approximation to the
%  coherent state \ket{\alpha} in the number basis.

% Ville Bergholm 2010


k = 0:n-1;
ket = (alpha.^k)./sqrt(factorial(k));
%ket = ket*exp(-abs(alpha)^2/2);
s = normalize(state(ket, n));

%s = exp(-abs(alpha)^2/2) * state(0, n).prop(expm(alpha*boson_ladder(n)'));
%s = state(0, n).prop(ho.displace(alpha, n));
