function [H, D] = fit(b, delta, T1, T2)
% FIT  Qubit-bath coupling that reproduces given decoherence times.
%
%   [H, D] = fit(b, delta, T1, T2)
%
%  Returns the single-qubit Hamiltonian H and the qubit-bath coupling operator D
%  that reproduce the decoherence times T1 and T2 (in TU)
%  for a single-qubit system coupled to the bath b.
%  delta is the energy splitting for the qubit (in units of \hbar/TU).
%
%  The bath object b is not modified in any way.

% Ville Bergholm 2009-2015


global qit;

if nargin ~= 4
  error('Required params: bath, delta, T1, T2.')
end

% T2 = 1/(0.5/T1 +1/T_dephase)
iTd = 1/T2 -0.5/T1; % inverse pure dephasing time
if iTd < 0
    error('Unphysical decoherence times!')
end

% match bath couplings to T1, T2
x = b.scale * delta/2;

switch b.type
  case 'ohmic'
    % Fitting an ohmic bath to a given set of decoherence times

    switch b.stat
      case 'boson'
        temp = coth(x) * x * b.cut_func(abs(delta));

        % coupling, ZX angle
        alpha = atan2(1, sqrt(T1 * iTd * temp));
        % dimensionless system-bath coupling factor squared
        c = iTd * b.scale/(4*pi*cos(alpha)^2);

        % noise coupling operator
        D = sqrt(c) * (cos(alpha)*qit.sz +sin(alpha)*qit.sx);

      case 'fermion'
        if iTd ~= 0
            error('For a fermionic bath we must have T2 = 2*T1')
        end
        % dimensionless system-bath coupling factor squared
        c = 1/(T1 * 2*pi * abs(delta) * b.cut_func(abs(delta)));
        D = sqrt(c) * qit.sx;
      otherwise
        error('Unknown statistics.')
    end

  otherwise
    error('Unsupported bath type.')
end

% qubit Hamiltonian
H = -delta/2 * qit.sz;
end
