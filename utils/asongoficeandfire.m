function y = asongoficeandfire(n)
% ASONGOFICEANDFIRE  Colormap with blues and reds. Wraps.
%  c = asongoficeandfire([n])
%
%  Returns an n * 3 matrix containing the colormap.
%  If n is not given, the colormap is the same length
%  as the current figure's colormap. 

% Ville Bergholm 2010


if (nargin < 1)
  n = size(get(gcf, 'colormap'), 1);
end

d = 3.1;
p = linspace(-1, 1, n);

x = p(find(p < 0));
y = [1-(1+x).^d; 0.5*(tanh(4*(-x-0.5))+1); (-x).^d];

x = p(find(p >= 0));
y = [y, [x.^d; 0.5*(tanh(4*(x-0.5))+1); 1-(1-x).^d]]';


return

figure
subplot(2,1,1)
plot(p, y')

subplot(2,1,2)
pcolor(kron(1:n, ones(10,1)))
shading flat
colormap(y);
end
