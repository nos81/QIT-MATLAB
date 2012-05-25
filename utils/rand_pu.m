function p = rand_pu(n)
% RAND_PU  Random n-partition of unity.
%  p = rand_pu(n)
%
%  p is random with respect to the order-n Dirichlet distribution Dir(\alpha)
%  with \alpha = (1, 1, ..., 1).

% Ville Bergholm 2008-2012


% Sample the Dirichlet distribution using n gamma-distributed
% random vars. The (shared) scale parameter of the gamma pdfs is irrelevant,
% and the shape parameters correspond to the Dirichlet \alpha params.
% Furthermore, Gamma(x; 1,1) = exp(-x), so

p = -log(rand(1, n)); % Gamma(1,1) -distributed
p = p / sum(p);  % Dir(1, 1, ..., 1) -distributed

% TODO this would be a simpler choice, but what's the exact distribution?
%x = sort(rand(n-1, 1));  % n-1 points in (0,1), sorted
%p = [x;1] -[0;x];        % n deltas between points = partition of unity
end
