function [g, s] = corr(b, x)
% BATH/CORR  Bath spectral correlation tensor.
%
%  [g, s] = corr(bath, dH)  % returns Gamma(omega0 * dH)/omega0, see below
%
%  Returns the bath spectral correlation tensor Gamma evaluated at omega0 * dH:
%
%    Gamma(omega0 * dH)/omega0 == 0.5*g +i*s

% Ville Bergholm 2009-2010


if (nargin ~= 2)
  error('Need a bath object and dH.')
end

tol = 1e-8;
max_w = 0.1; % TODO justify


% assume parameters are set and lookup table computed
%s = interp1(b.dH, b.s_table, x, 'linear', 0);

% binary search for the interval in which x falls
first = 1;
last = length(b.dH);
while (first+1 ~= last)
  pivot = round((first + last)/2);
  if (x < b.dH(pivot))
    last = pivot;
  else
    first = pivot;
  end
end
ee = b.dH([first last]); 
tt = b.gs_table(:, [first last]);
% now x is in [ee(1), ee(2))

gap = ee(2) - ee(1);
d1 = abs(x - ee(1));
d2 = abs(x - ee(2));

% close enough to either endpoint?
if (d1 <= tol)
  temp = b.gs_table(:, first);

elseif (d2 <= tol)
  temp = b.gs_table(:, last);

elseif (gap <= max_w + tol) % short enough gap to interpolate?
  temp = interpolate(ee, tt, x);

else % compute a new point
  if (gap <= 2*max_w)
    p = ee(1) +gap/2; % halfway
    if (x < p)
      idx = 2; % which ee p will replace
    else
      idx = 1;
    end
  elseif (d1 <= max_w)
    p = ee(1)+max_w;
    idx = 2;
  elseif (d2 <= max_w)
    p = ee(2)-max_w;
    idx = 1;
  else
    p = x;
    idx = 1;
  end

  % compute new g, S values at p and insert them into the table
  temp = [0; 0];
  temp(2) = S_func(b, p);

  if (abs(p) <= tol)
    temp(1) = b.g0; % limit at p == 0
  else
    temp(1) = b.g_func(p) .* b.cut_func(p);
  end

  b.dH = [b.dH(1:first), p, b.dH(last:end)];
  b.gs_table = [b.gs_table(:, 1:first), temp, b.gs_table(:, last:end)];

  % now interpolate the required value
  ee(idx) = p;
  tt(:, idx) = temp;
  temp = interpolate(ee, tt, x);
end

g = temp(1);
s = temp(2);

end

function y = interpolate(ee, tt, x)
% interp1 does way too many checks
  y = tt(:,1) + ((x - ee(1))/(ee(2) - ee(1)))*(tt(:,2) - tt(:,1));
end
