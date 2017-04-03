function [dH, A, Gamma] = components(H, D, baths)
% COMPONENTS  Components for the Lindblad equation.
%  [dH, A, Gamma] = components(H, D, B)

%  Builds the jump operators and computes the spectral correlation
%  tensor required for evaluating the Lindblad equation in the
%  weak-coupling approximation.
%
%  H is the system Hamiltonian, D a hermitian interaction operator
%  coupling the system to the bath B.
%
%  B can also be a cell vector of baths, in which case D has to be
%  a cell vector of the corresponding interaction operators.
%  The baths are assumed independent and uncorrelated,
%  i.e. \Gamma_{ab} = \delta_{ab} * \Gamma_{a}.

% Ville Bergholm 2017


if ~iscell(baths)
  baths = {baths}; % needs to be a cell array, even if it has just one element
end
n_baths = length(baths); % number of baths

% make sure the baths have the same TU!
temp = baths{1}.TU;
for k=2:n_baths
  if baths{k}.TU ~= temp
    error('All the baths must have the same time unit!')
  end
end

% transition freqs, jump ops
[dH, A] = rotating_frame(H, D);

if size(A,2) ~= n_baths
    error('Every bath must have a corresponding interaction operator and vice versa.')
end

n = length(dH);

% spectral correlation tensor evaluated at transition freqs
Gamma = zeros(n_baths, n, 2);
for k=1:n_baths
  b = baths{k};
  for j=1:n
      [g,s] = b.corr(dH(j));
      Gamma(k,j,1) = 0.5*g +1i*s;
      [g,s] = b.corr(-dH(j));
      Gamma(k,j,2) = 0.5*g +1i*s;
  end
end
