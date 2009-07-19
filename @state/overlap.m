function x = overlap(a, b)
% STATE/OVERLAP  Returns the overlap of two quantum states.
%  x = overlap(a, b)
%
%  For state vectors,   x = |<a|b>|^2.
%  For state operators, x = \trace(\rho_a * \rho_b)

% Ville Bergholm 2009


if (size(a.data, 2) == 1)
  if (size(b.data, 2) == 1)
    x = abs(a.data' * b.data)^2;
  else
    x = a.data' * b.data * a.data;
  end
else
  if (size(b.data, 2) == 1)
    x = b.data' * a.data * b.data;
  else
    x = trace(a.data * b.data);
  end
end
