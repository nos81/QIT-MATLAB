function c = canonical(U)
% CANONICAL  Canonical local invariants of a two-qubit gate.
%  c = canonical(U)
%
%  Returns a vector of three real canonical local invariants for the
%  U(4) matrix U, normalized to the range [0,1].

%! Childs et al., "Lower bounds on the complexity of simulating quantum gates", PRA 68, 052311 (2003).
% Ville Bergholm 2004-2010


global qit;

sigma = 1;
for k=1:2
  sigma = kron(sigma, qit.sy);
end

temp = U*sigma*U.'*sigma/sqrt(det(U));

lambda = eig(temp); %[exp(i*2*phi_1), etc]

% logarithm to the branch (-1/2, 3/2]
lambda = angle(lambda)/pi; % divide pi away
for k=1:4
  if (lambda(k) <= -1/2)
    lambda(k) = lambda(k) + 2;
  end
end
S = lambda/2;
S = sort(S, 1, 'descend');

n = round(sum(S)); % must be n

% take away extra translations-by-pi
S = S - cat(1, ones(n,1), zeros(4-n,1));
% put the elements in the correct order
S = circshift(S, -n);

M = [1 1 0; 1 0 1; 0 1 1]; % scaled by factor 2
c = (M*S(1:3)).';

% and into the Berkeley chamber using a translation and two Weyl reflections
if (c(3) < 0)
  c(1) = 1 - c(1);
  c(3) = -c(3);
end
c = mod(abs(c), 1);
