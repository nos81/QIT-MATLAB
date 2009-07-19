function quantum_channels(p)
% QUANTUM_CHANNELS  Visualization of simple one-qubit channels.
%
%  s = quantum_channels(p)
%
%  Visualizes the effect of different quantum channels on a qubit.

% Ville Bergholm 2009


fprintf('\n\n=== Quantum channels ===\n')

if (nargin < 1)
  p = 0.3
end

global qit;

I  = qit.I;
sx = qit.sx;
sy = qit.sy;
sz = qit.sz;

E_bitflip      = {sqrt(1-p)*I, sqrt(p)*sx};
E_phaseflip    = {sqrt(1-p)*I, sqrt(p)*sz};
E_bitphaseflip = {sqrt(1-p)*I, sqrt(p)*sy};
E_depolarize   = {sqrt(1-3*p/4)*I, sqrt(p)*sx/2, sqrt(p)*sy/2, sqrt(p)*sz/2};
%t = asin(sqrt(gamma))
t = pi/3;
E_amplitudedamp = {sqrt(p)*diag([1, cos(t)]), sqrt(p)*[0 sin(t); 0 0], sqrt(1-p)*diag([cos(t), 1]), sqrt(1-p)*[0 0; sin(t) 0]};

[X,Y,Z] = sphere(20);
S(1,:,:) = X;
S(2,:,:) = Y;
S(3,:,:) = Z;

%present(S, E_bitflip, 'Bit flip channel');
%present(S, E_depolarize, 'Depolarizing channel');
present(S, E_amplitudedamp, 'Amplitude damping channel');
end



function res = present(S, E, T)

s = size(S);
res = [];

for a=1:s(2)
  for b=1:s(3)
    temp = bloch_state(S(:,a,b));
    res(:,a,b) = bloch_vector(kraus_propagate(temp, E));
  end
end

figure
plots.bloch_sphere();
surf(squeeze(res(1,:,:)), squeeze(res(2,:,:)), squeeze(res(3,:,:)));
shading flat
alpha(0.2)
title(T);
end
