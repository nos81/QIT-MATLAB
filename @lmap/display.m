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

if (n_ind == 1 || (n_ind == 2 && (isequal(s.dim{1}, 1) || isequal(s.dim{2}, 1))))
  % vector or a map with a singleton domain or codomain dimension

  % ket or bra?
  if (n_ind == 1 || isequal(s.dim{2}, 1))
    dim = s.dim{1};
    is_ket = true;
  else
    dim = s.dim{2};
    is_ket = false;
  end

  if (length(s.data) > 128)
    % sanity check, do not display lmaps with hundreds of terms
    out = [out, ' (long)'];
  else

  n = length(dim);

  for ind = 1:prod(dim)
    temp = s.data(ind);
    if (abs(temp) >= qit.tol)
      if (abs(imag(temp)) <= qit.tol)
        out = [out, sprintf(' %+2g', real(temp))];
      elseif (abs(real(temp)) <= qit.tol)
        out = [out, sprintf(' %+2gi', imag(temp))];
      else
        out = [out, sprintf(' +(%s)', num2str(temp))];
      end
      
      d = fliplr(dim); % big-endian convention makes this way more complicated than it should be
      carry = ind;
      for k = 1:n % start from least significant digit
        %ind2sub with 2 output parms uses up the first dim given
        [ket(k), carry] = ind2sub(d(k:n), carry);
      end
      ket = fliplr(ket); % big-endian again

      % ket or bra?
      if (is_ket)
        out = [out, ' |', char(ket - 1 + '0'), '>'];
      else
        out = [out, ' <', char(ket - 1 + '0'), '|'];
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
    out = [out, '[ ', sprintf('%d ', s.dim{k}), ']  '];
    if (k < n_ind)
      out = cat(2, out, '<-  ');
    end
  end
  disp(out);
else
  disp(' ') % empty line
end
