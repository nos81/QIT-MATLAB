function L = superop(H, D, baths)
% SUPEROP  Liouvillian superoperator for a Born-Markov system.
%
%  L = superop(H, D, B)
%
%  Builds the Liouvillian superoperator L corresponding to a
%  system Hamiltonian H and a (Hermitian) interaction operator D
%  coupling the system to bath B in the weak coupling limit.
%
%  Returns L * TU, which includes the system Hamiltonian, the Lamb shift,
%  and the Lindblad dissipator.
%
%  B can also be a cell vector of baths, in which case D has to be
%  a cell vector of the corresponding interaction operators.
%  The baths are assumed independent and uncorrelated,
%  i.e. \Gamma_{ab} = \delta_{ab} * \Gamma_{a}.

% Ville Bergholm 2009-2017


[dH, X, Gamma] = markov.components(H, D, baths);
n_baths = size(X,2);

iH_LS = 0;
acomm = 0;
diss = 0;

for n=1:n_baths
  A = X(:,n); % ops for the n'th bath

  % we build the Liouvillian in a funny order to be a bit more efficient
  for k=1:length(dH)
    % first the positive energy shift
    G = Gamma(n, k, 1);
    g = 2*real(G);
    temp = A{k}' * A{k};

    iH_LS = iH_LS +1i * imag(G) * temp;
    acomm = acomm -0.5 * g * temp;
    diss = diss +lrmul(g * A{k}, A{k}');

    if dH(k) == 0
        % no negative shift
        continue
    end

    % now the corresponding negative energy shift
    G = Gamma(n, k, 2);
    g = 2*real(G);
    temp = A{k} * A{k}'; % note the difference here, A(-omega) = A'(omega)

    iH_LS = iH_LS +1i * imag(G) * temp;
    acomm = acomm -0.5 * g * temp;
    diss = diss +lrmul(g * A{k}', A{k}); % here too
  end
end

iH_LS = iH_LS +1i*H; % include the system Hamiltonian

L = lmul(acomm -iH_LS) +rmul(acomm +iH_LS) +diss;
