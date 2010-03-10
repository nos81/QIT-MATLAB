function display(s, short)
% STATE/DISPLAY  Print the state in a neat format.

% Ville Bergholm 2008-2010


global qit;

if (nargin < 2)
  short = false;
end

out = [inputname(1) ' ='];

if (~short)
  % long format
  disp(out);
  disp('state:');
  out = '  ';
end

if (size(s.data, 2) == 1)
  % state vector (pure state)

  if (size(s.data, 1) > 128)
    % sanity check, do not display states with hundreds of terms
    out = cat(2, out, ' (long)');
  else

  n = length(s.dim);

  for ind = 1:prod(s.dim)
    temp = s.data(ind);
    if (abs(temp) >= qit.tol)
      if (abs(imag(temp)) <= qit.tol)
        out = cat(2, out, sprintf(' %+2g', real(temp)));
      elseif (abs(real(temp)) <= qit.tol)
        out = cat(2, out, sprintf(' %+2gi', imag(temp)));
      else
        out = cat(2, out, sprintf(' +(%s)', num2str(temp)));
      end
      
      d = fliplr(s.dim); % big-endian convention makes this way more complicated than it should be
      carry = ind;
      for k = 1:n % start from least significant digit
        %ind2sub with 2 output parms uses up the first dim given
        [ket(k), carry] = ind2sub(d(k:n), carry);
      end
      ket = fliplr(ket); % big-endian again

      out = cat(2, out, ' |');
      out = cat(2, out, char(ket - 1 + '0'));
      out = cat(2, out, '>');
    end
  end
  end
  disp(out);
else
  % state operator
  disp(s.data);
end

if (~short)
  disp('dim:');
  disp(s.dim);
else
  disp(' ') % empty line
end
