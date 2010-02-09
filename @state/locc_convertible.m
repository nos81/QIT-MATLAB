function res = locc_convertible(s, t, sys)
% STATE/LOCC_CONVERTIBLE  LOCC convertibility of states.
%  res = locc_convertible(s, t, sys)
%
%  For bipartite pure states s and t, returns true if s can be converted to t
%  using local operations and classical communication (LOCC).
%  sys is a vector of subsystems defining the partition.

% Ville Bergholm 2010
%! M.A. Nielsen, I.L. Chuang, "Quantum Computation and Quantum Information" (2000), chapter 12.5.1


if (size(s.data, 2) ~= 1)
  s = to_ket(s);
  %error('Not implemented for nonpure states.')
end

if (size(t.data, 2) ~= 1)
  t = to_ket(t);
  %error('Not implemented for nonpure states.')
end

if (length(s.dim) ~= length(t.dim) || any(s.dim ~= t.dim))
  error('States must have equal dimensions.')
end

s = ptrace(s, sys);
t = ptrace(t, sys);

res = majorize(eig(s.data), eig(t.data));