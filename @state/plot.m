function h = plot(s)
% PLOT  State tomography plot.
%  h = plot(s)
%
%  Plots the probabilities of finding a system in the state s
%  in the different computational basis states upon measurement.
%  Relative phases are represented by the colors of the bars.
%
%  If s is a nonpure state, also plots the coherences.

% Ville Bergholm 2009-2010


% prepare labels
dim = dims(s);
n = length(dim);

m = min(n, 3); % at most three symbols
d = dim(1:m);
nd = prod(d);

ntot = prod(dim);
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
Ncol = 127; % color resolution (odd to fix zero phase at the center of a color index)
colormap(circshift(hsv(Ncol), floor(Ncol/6))); % the hsv map wraps (like phase)

if is_ket(s)
  % state vector
  s = fix_phase(s);

  h = bar(0:N-1, prob(s));
  xlabel('Basis state');
  ylabel('Probability');
  set(gca,'XTick', ticks)
  set(gca,'XTickLabel', ticklabels);
  axis tight;

  % color bars using phase data
  ch = get(h,'Children');
  fvd = get(ch,'Faces');
  fvcd = get(ch,'FaceVertexCData');

  c = phases(s.data);
  for b = 1:N
    fvcd(fvd(b,:)) = c(b); % all four vertices of a bar have same color
  end
  set(ch,'FaceVertexCData',fvcd);
  set(ch,'EdgeColor','k');

else
  % state op
  h = bar3(abs(s.data));
  xlabel('Col state');
  ylabel('Row state');
  zlabel('|\rho|');
  set(gca,'XTick', ticks+1);
  set(gca,'XTickLabel', ticklabels);
  set(gca,'YTick', ticks+1);
  set(gca,'YTickLabel', ticklabels);
  axis tight;
  %alpha(0.8);

  % color bars using phase data
  c = phases(s.data);
  for m = 1:length(h)
    % get color data
    cdata = get(h(m), 'Cdata'); % [one row of 3d bars * six faces, four vertices per face]
    for k = 1:(size(cdata, 1)/6)
      j = 6*k;
      cdata(j-5:j, :) = c(k,m); % all faces are the same color
    end
    set(h(m), 'Cdata', cdata);
  end

end

set(gca, 'CLim', [0 1]); % color limits

hcb = colorbar('YTick', linspace(0, 1, 5));
set(hcb, 'YTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
end


function p = phases(A)
  p = 0.5*((angle(A)/pi)+1); % normalized to (0,1]
end


function c = colors(p, n)
  % p in [0,1) is uniformly mapped to {1, 2, ..., n}
  %c = 1 + floor(n * p);

  % p in (0,1] is uniformly mapped to {1, 2, ..., n}
  c = ceil(n * p);
end
