function display(s)
% STATE/DISPLAY  Display the state.

% Ville Bergholm 2008-2009


global qit;

if (size(s.data, 2) == 1)
  % state vector (pure state)

  out = 'state:  ';
  n = length(s.dim);

  for ind = 1:prod(s.dim)
    temp = s.data(ind);
    if (abs(temp) >= qit.tol)
      if (abs(imag(temp)) <= qit.tol)
        out = strcat(out, sprintf(' %+2g', real(temp)));
      elseif (abs(real(temp)) <= qit.tol)
        out = strcat(out, sprintf(' %+2gi', imag(temp)));
      else
        out = strcat(out, sprintf(' +(%s)', num2str(temp)));
      end
      
      d = fliplr(s.dim); % big-endian convention makes this way more complicated than it should be
      carry = ind;
      for k = 1:n % start from least significant digit
        %ind2sub with 2 output parms uses up the first dim given
        [ket(k), carry] = ind2sub(d(k:n), carry);
      end
      ket = fliplr(ket); % big-endian again

      out = strcat(out, ' |');
      out = strcat(out, char(ket - 1 + '0'));
      out = strcat(out, '>');
    end
  end
  disp(out);
else
  % state operator
  disp('state:');
  disp(s.data);
end

disp('dim:');
disp(s.dim);
