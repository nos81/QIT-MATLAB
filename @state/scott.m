function Q = scott(s, m)
% SCOTT  Scott's average bipartite entanglement measure.
%  Q = scott(s, m)
%
%  Returns the vector Q containing the terms of the Scott entanglement measure
%  of system s for partition size m.
%
%  When m = 1 this is coincides with the Meyer-Wallach entanglement measure.

%! P.J. Love et al., "A characterization of global entanglement", arXiv:quant-ph/0602143 (2006).
%! A.J. Scott, "Multipartite entanglement, quantum-error-correcting codes, and entangling power of quantum evolutions", PRA 69, 052330 (2004).
%! D.A. Meyer and N.R. Wallach, "Global entanglement in multiparticle systems", J. Math. Phys. 43, 4273 (2002).
%
% Jacob D. Biamonte 2008
% Ville Bergholm 2008-2010


dim = dims(s);
n = length(s.dim);    % number of subsystems

S = nchoosek(1:n, m); % all m-combinations of n subsystems
n_c = size(S, 1); % number of distinct combinations

D = min(dim); % FIXME correct for arbitrary combinations of qudits??
C = (D^m/(D^m-1))/nchoosek(n, m); % normalization constant

Q = [];
for k = 1:n_c
  temp = ptrace(s, setdiff(1:n, S(k, :))); % trace over everything except S_k
  % NOTE: For pure states, tr(\rho_S^2) == tr(\rho_{\bar{S}}^2),
  % so for them we could just use ptrace(s, S(k, :)) here.
  temp = temp.data;
  Q(k) = C*(1 - trace(temp^2)); 
end
