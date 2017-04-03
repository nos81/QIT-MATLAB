function [LL, H_LS] = lindblad_ops(H, D, baths)
% LINDBLAD_OPS  Lindblad operators for a Born-Markov master equation.
%  [LL, H_LS] = lindblad_ops(H, D, B)
%
%  Builds the Lindblad operators corresponding to a
%  system Hamiltonian H and a (Hermitian) interaction operator D
%  coupling the system to bath B in the weak coupling limit.
%
%  Returns LL == {L_i * sqrt(TU)}_i and H_LS * TU / \hbar,
%  where L_i are the Lindblad operators and H_LS is the Lamb shift.
%
%  B can also be a cell vector of baths, in which case D has to be
%  a cell vector of the corresponding interaction operators.

% Ville Bergholm 2009-2017


global qit


[dH, X, Gamma] = markov.components(H, D, baths);
n_baths = size(X,2);

H_LS = 0;

ind = 1;
for n=1:n_baths
  A = X(:,n); % ops for the n'th bath

  for k=1:length(dH)
    % first the positive energy shift
    G = Gamma(n, k, 1);
    g = 2*real(G);

    % is the dissipation significant?
    temp = sqrt(g) * A{k};
    if norm(temp) >= qit.tol
        LL{ind} = temp;
        ind = ind+1;
    end
    % contribution to Lamb shift
    H_LS = H_LS +imag(G) * A{k}' * A{k};

    if dH(k) == 0
        % no negative shift
        continue
    end

    % now the corresponding negative energy shift
    G = Gamma(n, k, 2);
    g = 2*real(G);

    temp = sqrt(g) * A{k}';   % note the difference here, A(-omega) = A'(omega)
    if norm(temp) >= qit.tol
        LL{ind} = temp;
        ind = ind+1;
    end
    H_LS = H_LS +imag(G) * A{k} * A{k}'; % here too
  end
end

% TODO ops for different baths can be combined into a single basis,
% N^2-1 ops max in total
