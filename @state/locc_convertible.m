function res = locc_convertible(s, t, sys)
% LOCC_CONVERTIBLE  LOCC convertibility of states.
%  res = locc_convertible(s, t, sys)
%
%  For bipartite pure states s and t, returns true if s can be converted to t
%  using local operations and classical communication (LOCC).
%  sys is a vector of subsystems defining the partition.

% Ville Bergholm 2010-2011
%! M.A. Nielsen, I.L. Chuang, "Quantum Computation and Quantum Information" (2000), chapter 12.5.1


s = to_ket(s);
t = to_ket(t);
%error('Not implemented for nonpure states.')

if (~equal_dims(s, t))
  error('States must have equal dimensions.')
end

s = ptrace(s, sys);
t = ptrace(t, sys);

res = majorize(eig(s.data), eig(t.data));
