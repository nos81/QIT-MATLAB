function display(s, short)
% DISPLAY  Display the lmap in a neat format.

% Ville Bergholm 2008-2010


global qit;

if (nargin < 2)
  short = false;
end

out = [inputname(1) ' ='];

if (~short)
  % long format
  disp(out);
  out = '  ';
end

% number of indices
n_ind = order(s);

if (n_ind == 1 || (n_ind == 2 && isequal(s.dim{1}, 1)))
  % vector

  % ket or bra?
  q = 2;
  if (n_ind == 1)
    q = 1;
  end

  if (length(s.data) > 128)
    % sanity check, do not display lmaps with hundreds of terms
    out = cat(2, out, ' (long)');
  else

  n = length(s.dim{q});

  for ind = 1:prod(s.dim{q})
    temp = s.data(ind);
    if (abs(temp) >= qit.tol)
      if (abs(imag(temp)) <= qit.tol)
        out = cat(2, out, sprintf(' %+2g', real(temp)));
      elseif (abs(real(temp)) <= qit.tol)
        out = cat(2, out, sprintf(' %+2gi', imag(temp)));
      else
        out = cat(2, out, sprintf(' +(%s)', num2str(temp)));
      end
      
      d = fliplr(s.dim{q}); % big-endian convention makes this way more complicated than it should be
      carry = ind;
      for k = 1:n % start from least significant digit
        %ind2sub with 2 output parms uses up the first dim given
        [ket(k), carry] = ind2sub(d(k:n), carry);
      end
      ket = fliplr(ket); % big-endian again

      % ket or bra?
      if (q == 1)
        out = cat(2, out, ' |', char(ket - 1 + '0'), '>');
      else
        out = cat(2, out, ' <', char(ket - 1 + '0'), '|');
      end
    end
  end
  end
  disp(out);
else
  % matrix or higher-order tensor
  disp(s.data);
end


if (~short)
  disp('dim:');
  out = '  ';
  for k = 1:n_ind
    out = cat(2, out, '[ ', sprintf('%d ', s.dim{k}), ']  ');
    if (k < n_ind)
      out = cat(2, out, '<-  ');
    end
  end
  disp(out);
else
  disp(' ') % empty line
end
