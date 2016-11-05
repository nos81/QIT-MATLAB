function display(s, no_dims, labels)
% DISPLAY  Display the lmap in a neat format.

% Ville Bergholm 2008-2016


global qit;

if nargin < 2
  no_dims = false;
end

out = [inputname(1) ' ='];

if ~no_dims
  % long format
  disp(out);
  out = '  ';
end

% number of indices
n_ind = order(s);

if n_ind == 1 || (n_ind == 2 && (isequal(s.dim{1}, 1) || isequal(s.dim{2}, 1)))
  % vector, or a map with a singleton domain or codomain dimension

  % ket or bra?
  if n_ind == 1 || isequal(s.dim{2}, 1)
    dim = s.dim{1};
    is_ket = true;
  else
    dim = s.dim{2};
    is_ket = false;
  end

  % TEST define state labels
  if nargin < 3
    temp = max(dim);
    labels = char((1:temp) -1 +'0');
  end
  
  D = prod(dim);
  if nnz(s.data) > 128
    % sanity check, do not display lmaps with hundreds of terms
    out = [out, ' (long)'];
  
  elseif D == 1
    % scalar, just print the value
    out = [out, sprintf(' %g', full(s.data))];
    no_dims = true;
    
  else
    % vector, print ket (or bra) symbols
    n = length(dim);

    for ind = 1:D
      temp = full(s.data(ind)); % sprintf can't handle sparse arrays, not even scalars
      if abs(temp) >= qit.tol
        if abs(imag(temp)) <= qit.tol
          out = [out, sprintf(' %+2g', real(temp))];
        elseif abs(real(temp)) <= qit.tol
          out = [out, sprintf(' %+2gi', imag(temp))];
        else
          out = [out, sprintf(' +(%s)', num2str(temp))];
        end

        ket = unravel_index(ind, dim);
        symbol = labels(ket);

        % ket or bra?
        if is_ket
          out = [out, ' |', symbol, '>'];
        else
          out = [out, ' <', symbol, '|'];
        end
      end
    end
  end
  disp(out);
else
  % matrix or higher-order tensor
  disp(s.data);
end


if ~no_dims
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
