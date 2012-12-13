function [block] = tdmrg(block, H_func, t, steps, m)
% TDMRG  Time-dependent DMRG.
%  B = tdmrg(B, H_func, t, m)
%
%  Propagates the DMRG block structure B in time using Hamiltonian
%  given by H_func.
%
%  [H, dim] = H_func([first, last]) should return the partial
%  Hamiltonian and dimension vector for the sites first:last.
%  [S1, S2, h1, h2] = H_func(s) should return the coupling
%  and local operators for sites s and s+1. 
%
%  t is the total propagation time.
%  At each truncation step, at most m states are kept (possibly
%  fewer near the ends of the chain).

%! G. De Chiara et al., "Density Matrix Renormalization Group for Dummies", JCTN 5, 1277 (2008), doi:10.1166/jctn.2008.011

% Ville Bergholm 2010


if (nargin < 5)
  m = 15
  if (nargin < 4)
    error('block, H_func, t and steps are required.');
  end
end

n = length(block)+1

dt = t/steps;

for k=1:steps
  % do a back-and-forth sweep
  q = 1;
  sweep_right = true;
  while (1)
    % after this iteration, block q will be replaced
    
    % local Hamiltonian terms (**)
    [block{q}.S1, block{q}.S2, block{q}.h1, block{q}.h2] = H_func(q);

    [H, dim] = starstar(block{q});
    
    if (q == 1)
      % symmetric sum, h1/2, h2/2, extra halves at ends
      H = H +kron(block{q}.h1/2, speye(dim(2)));
    end
    
    if (q == n-1)
      % last block, reverse direction
      sweep_right = false;
      
      %H = H +kron(speye(dim(1)), block{q}.h2); % add h_n term(this is left-right unsymmetric!)

      % symmetric sum, h1/2, h2/2, extra halves at ends
      H = H +kron(speye(dim(1)), block{q}.h2/2);

      ddt = dt; % full step (since we do it just once)
    else
      ddt = dt/2; % half step
    end
    H = 0.5*(H+H'); % eliminate rounding errors

    % propagate state
    dim = block{q}.state.dim;
    block{q}.state = u_propagate(block{q}.state, mkron(speye(dim(1)), expm(full((-1i*ddt)*H)), speye(dim(4))));

    % truncate
    block{q} = truncate(block{q}, m);

    % project state to next block using the state prediction algorithm
    if (sweep_right)
      block{q+1}.state = guess_stateR(block{q}, block{q+1});
      q = q+1;
    else
      if (q == 1)
        break; % done
      end
      block{q-1}.state = guess_stateL(block{q}, block{q-1});
      q = q-1;
    end
  end
end

end


% White's state prediction
function s = guess_stateR(B, C)
% make a guess towards the right
  dim = B.state.dim;
  O1 = B.OL';
  O2 = C.OR;
  c = length(C.h2);
  dd = [size(O1,1), length(C.h1), c, size(O2,1)/c];
  s = state(kron(speye(size(O1, 1)*dim(3)), O2)*(kron(O1, speye(prod(dim(3:4))))*B.state.data), dd);
end
      
function s = guess_stateL(B, A)
  dim = B.state.dim;
  O1 = B.OR';
  O2 = A.OL;
  a = length(A.h1);
  dd = [size(O2,1)/a, a, length(A.h2), size(O1,1)];
  s = state(kron(O2, speye(dim(2)*size(O1, 1)))*(kron(speye(prod(dim(1:2))), O1)*B.state.data), dd);
end

function C = coupling(S1, S2)
% C = \sum_k S1_k * S2_k
  C = sparse(0);
  for k=1:length(S1)
    C = C +kron(S1{k}, S2{k});
  end
end


function [H, dim] = starstar(B)
% ** block
  dim = [length(B.h1), length(B.h2)];

  % TEST only include h1 here, h2 only in last cell
  %H = kron(B.h1, speye(dim(2))) +coupling(B.S1, B.S2);% +kron(speye(dim(1)), B.h2);

  % symmetric version, include only half of local terms at each block
  H = kron(B.h1/2, speye(dim(2))) +coupling(B.S1, B.S2) +kron(speye(dim(1)), B.h2/2);
end


function B = truncate(B, m)
  % Schmidt decompose block state (4 components) into two halves
  [d, u, v] = schmidt(B.state, [1 2]);

  keep = min(m, length(d)); % how many states to keep?
    
  B.P = 1 - sum(d(1:keep).^2); % truncation error
  % keep most influential states
  B.OL = u(:, 1:keep);
  B.OR = v(:, 1:keep);    
end
