function G = gellmann(n)
% GELLMANN  Gell-Mann matrices of dimension n.
%  G = gellmann(n)
%
%  Returns the n^2-1 (traceless, Hermitian) Gell-Mann matrices of dimension n,
%  normalized such that trace(G_i' * G_j) = \delta_ij.

% Ville Bergholm 2006-2009


global qit;

if (n <= 0)
  error('Dimension must be greater than one.')
end

% check cache first
if (length(qit.gellmann) >= n && length(qit.gellmann{n}) > 0)
  G = qit.gellmann{n};
  return;
end


count = 1;
d = zeros(1,n);
d(1) = 1;
for k=2:n
  for j=1:k-1
    temp = zeros(n);
    temp(k,j) = 1;
    temp(j,k) = 1;
    G{count} = temp/sqrt(2);
  
    temp(k,j) = i;
    temp(j,k) = -i;
    G{count+1} = temp/sqrt(2);

    count = count+2;
  end
  d(k) = -sum(d);
  G{count} = diag(d)/sqrt(sum(d.^2));
  count = count+1;
  d(k) = 1; 
end

% store them in the cache
qit.gellmann{n} = G;
