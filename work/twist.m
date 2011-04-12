function T = twist(s, sys, sig)
% TWIST  Twist invariants of an n-qubit state.
%  T = twist(s, sys, sig)
%
% Given an n-qubit state s, returns the twist T in the loop consisting of the subsystems
% in the vector sys. If sys is not given, it is assumed to be 1:n.
%
% The vector sig contains the loop signature used. Default is all ones.

%! M. S. Williamson et al., "Quantum correlations can be twisted", arXiv:1102.5609 (2011). 
% Ville Bergholm 2011


n = length(s.dim);

if (nargin < 2)
  % default: loop over all systems
  sys = 1:n;
end

m = length(sys); % loop length

if (nargin < 3)
  % which signatures should we use?
  sig = 1*ones(m,1);
end

eta = diag([1, -1, -1, -1]);
% allowed signature transforms (\in SO^+(1,3))
Q = {diag([1 1 1 1]), diag([1 -1 -1 1]), diag([1 -1 1 -1]), diag([1 1 -1 -1])};
swap = gate.swap(2,2); % 2q swap gate


T = eye(4);
err = false;

for k = 1:m
  a = sys(k);
  b = sys(mod(k,m) + 1);
  ignore = setdiff(1:n, [a b]); % trace out everything except a and b
  temp = ptrace(s, ignore);
  % always make sure system b (the one the local op acts on) is on the right
  if (b < a)
    temp = u_propagate(temp, swap);
  end
  [V, Sigma, W, S, errflag] = lorentz_svd(temp);
  err = err || errflag;

  % symmetrize on the right side: S = \tilde{S} * \Lambda
  Lambda = eta * V * Q{sig(k)} * eta * W.';
  
  T = T*Lambda; % multiply from the right so T is SLOCC-invariant
  %T = Lambda*T; % wrong order
end

T = 0.25*trace(T);
if (err)
    T = -1;
end
end
