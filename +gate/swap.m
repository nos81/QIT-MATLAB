function U = swap(d1, d2)
% SWAP  Swap gate.
%  U = swap(d1, d2)
%
%  Returns the swap gate which swaps the order of two subsystems with dimensions [d1 d2].
%  S: A_1 \otimes A_2 \to A_2 \otimes A_1,
%  S(v_1 \otimes v_2) = v_2 \otimes v_1  for all v_1 \in A_1, v_2 \in A_2.

% Ville Bergholm 2010


temp = d1*d2;
U = sparse(temp, temp);
for x = 0:d1-1
  for y = 0:d2-1
    U(d1*y+x+1, d2*x+y+1) = 1;
  end
end

U = lmap(U, {[d2 d1], [d1 d2]});
