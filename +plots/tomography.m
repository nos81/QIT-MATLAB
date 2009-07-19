function tomography(s)
% PLOT/TOMOGRAPHY  State tomography plot.
%  tomography(s)
%
%  Plots the probabilities of finding a system in the state s
%  in different computational basis states upon measurement.
%  If s is a nonpure state, also plots the coherences.

% Ville Bergholm 2009


% prepare labels
n = length(s.dim);
m = min(n, 3); % at most three symbols
d = s.dim; d = d(1:m); % stupid HACK
nd = prod(d);

ntot = prod(s.dim);
skip = ntot/nd; % only every skip'th state gets a label to avoid clutter
ticks = 0:skip:ntot-1;

ket = zeros(1,n);
ticklabels = {};
for k=0:nd-1
  ticklabels{end+1} = char(ket + '0');

  for b = m:-1:1 % start from least significant digit
    ket(b) = ket(b)+1;
    if (ket(b) < d(b))
      break;
    end
    ket(b) = 0;
  end
end


N = size(s.data, 1);

if (size(s.data, 2) == 1)
  % state vector
  bar(0:N-1, prob(s));
  xlabel('Basis state');
  ylabel('Probability');
  set(gca,'XTick', ticks)
  set(gca,'XTickLabel', ticklabels)
else
  bar3(abs(s.data))
  xlabel('Basis state');
  ylabel('Basis state');
  zlabel('|\rho|');
  %title('Absolute value of \rho')
  set(gca,'XTick', ticks+1)
  set(gca,'XTickLabel', ticklabels)
  set(gca,'YTick', ticks+1)
  set(gca,'YTickLabel', ticklabels)
  %alpha(0.6)
  % subplot(3,4,5)
  %bar3(angle(rho(:,:,p)))
  %axis([0 5 0 5 -pi pi])
  %title('Phase of \rho')
end
