function display(s)
% STATE/DISPLAY  Display the state.

% Ville Bergholm 2008


if (size(s.data, 2) == 1)
  % pure state (state vector)

  out = 'state:  ';
  n = length(s.dim);

  for ind = 1:prod(s.dim)
    if (s.data(ind) ~= 0)
      out = strcat(out, sprintf(' %+2g', s.data(ind)));

      d = fliplr(s.dim); % FUCK, big-endian convention makes this way more complicated than it should be
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
  disp('state:');
  disp(s.data);
end

disp('dim:');
disp(s.dim);
