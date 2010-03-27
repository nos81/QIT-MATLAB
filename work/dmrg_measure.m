function res = dmrg_measure(block, oplist)
% DMRG_MEASURE  DMRG measurement.
%  res = dmrg_measure(B, oplist)
%
%  Returns the approximate expectation values of the operators
%  defined in oplist (see op_list documentation for the format
%  definition) in the cell vector res.
%
%  B is a DMRG block struct which contains the (converged)
%  block states, obtained from the function dmrg_finite. 

% Ville Bergholm 2010


% TODO for the moment we only measure local or nn-observables

n = length(block)+1;

% loop over all operator definitions in the list
for k=1:length(oplist)
  op = oplist{k};
  no = size(op, 1); % number of constituent ops

  q = op{1,2}; % site of first op
  if (no == 1)
      if (q == n)
        q = n-1;
        dim = block{q}.state.dim;
        X = mkron(speye(prod(dim(1:2))), op{1,1}, speye(dim(4)));
      else
        dim = block{q}.state.dim;
        X = mkron(speye(dim(1)), op{1,1}, speye(prod(dim(3:4))));
      end
      
  elseif (no == 2)
      if (op{2,2} ~= q+1)
        error('not nn')
      end

      dim = block{q}.state.dim;
      X = mkron(speye(dim(1)), op{1,1}, op{2,1}, speye(dim(4)));
    
    else
      error('TODO FIXME');
  end
  
  res{k} = ev(block{q}.state, X);
end
