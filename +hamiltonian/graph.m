function [H, dim] = graph(G)
% HAMILTONIAN/GRAPH  Two-local Hamiltonians from connection graphs.
%  H = graph(G)         % just the Hamiltonian
%  [H, dim] = graph(G)  % Hamiltonian and dimension vector
%
%  Builds a Hamiltonian H corresponding to the graph G.
%  G is given in the form of an upper triangular connection
%  matrix array, where each diagonal element (k,k) is a
%  list of local terms operating on the k:th subsystem,
%  and each nondiagonal element (j,k) a list of pairs, where
%  each pair {A,B} corresponds to the term A_j*B_k.
%  Hence H consists of at most two-local terms.
%
%  Example: The connection graph
%    G = {{sz}, {},   {{sx,sx}, {sy,sy}, {sz,sz}}; ...
%         {},   {sz}, {{A,B}}; ...
%         {},   {},   {2*sz}}
%  corresponds to the Hamiltonian
%    H = sz_1 +sz_2 +2*sz_3 +sx_1*sx_3 +sy_1*sy_3 +sz_1*sz_3 +A_2*B_3.
%
%  The output parameter dim contains the dimensionality of the system
%  H is operating on, inferred from G.

% Ville Bergholm 2009


n = length(G);
H = 0;

for a=1:n
  % 1-local terms
  % allowing several terms here is useless (you could just sum them up), but at least we're consistent
  q = length(G{a,a});
  if (q > 0)
    temp = 0;
    for c=1:q
      temp = temp +G{a,a}{c};
    end
    H = H +mkron(eye(2^(a-1)), temp, eye(2^(n-a)));
  end
    
  for b=a+1:n
    % 2-local terms
    q = length(G{a,b});
    if (q > 0)
      temp = 0;
      for c=1:q
        temp = temp +mkron(G{a,b}{c}{1}, eye(2^(b-a-1)), G{a,b}{c}{2});
      end
      H = H +mkron(eye(2^(a-1)), temp, eye(2^(n-b)));
    end
  end
end

if (nargout == 2)
  dim = zeros(1,n);
  for a=1:n
    dim(a) = length(G{a,a});
  end
end
