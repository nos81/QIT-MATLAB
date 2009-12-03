function x = overlap(a, b)
% STATE/OVERLAP  Overlap of two states.
%  x = overlap(a, b)
%
%  For state vectors,   x = |<a|b>|.
%  For state operators, x = \sqrt{\trace(\rho_a * \rho_b)}

% Ville Bergholm 2009


if (size(a.data, 2) == 1)
  if (size(b.data, 2) == 1)
    x = abs(a.data' * b.data);
  else
    x = sqrt(a.data' * b.data * a.data);
  end
else
  if (size(b.data, 2) == 1)
    x = sqrt(b.data' * a.data * b.data);
  else
    x = sqrt(trace(a.data * b.data));
  end
end
