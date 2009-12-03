function x = ev(s, A)
% EV  Expectation value of an observable in the state.
%  x = ev(s, A)
%
%  Returns the expectation value of observable A in the state s.
%  A has to be Hermitian.

% Ville Bergholm 2008


if (size(s.data, 2) == 1)
  % state vector
  x = s.data' * A * s.data;
else
  % state operator
  x = trace(A*s.data);
end

x = real(x); % Hermitian observable
