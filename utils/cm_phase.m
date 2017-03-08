function y = cm_phase(n)
% CM_PHASE  Colormap for denoting complex phases.
%  c = cm_phase([n])
%
%  Returns a colormap that represents a closed path (more or less a circle) in RGB space.
%
%  Returns an n * 3 matrix containing the colormap.
%  If n is not given, the colormap is the same length
%  as the current figure's colormap.

% Ville Bergholm 2016-2017


if nargin < 1
  n = size(get(gcf, 'colormap'), 1);
end

% the hsv map wraps so we could just use that
%y = circshift(hsv(n), floor(n/6));
%return

p = linspace(-pi, pi, n+1);
p = p(2:end);

r = 0.5 * [1; 1; 1];  % RGB cube center
v = [-1; -1; 2]/4;
% v is orthogonal to r, r+v touches the mostly yellow surface
a = [1;1;1]/sqrt(3); % rotation axis
for k=1:n
    y(k, :) = r +R_SO3(a * p(k)) * v;
end
% clamp to the RGB cube
y = min(y, 1);
y = max(y, 0);
return

figure
subplot(2,1,1)
plot(p, y')

subplot(2,1,2)
pcolor(kron(1:n, ones(10,1)))
shading flat
colormap(y);
end
