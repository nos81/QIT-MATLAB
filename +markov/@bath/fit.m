function [H, D] = fit(b, delta, T1, T2)
% FIT  Qubit-bath coupling that reproduces given decoherence times.
%
%   [H, D] = fit(b, delta, T1, T2)
%
%  Returns the qubit Hamiltonian H and the qubit-bath coupling operator D
%  that reproduce the decoherence times T1 and T2 (in units of 1/omega0)
%  for a single-qubit system coupled to the bath b.
%  delta is the energy splitting for the qubit (in units of hbar*omega0).
%
%  The bath object b is not modified in any way.

% Ville Bergholm 2009-2010


global qit;

if (nargin ~= 4)
  error('Required params: bath, delta, T1, T2.')
end

switch b.type
  case 'ohmic'
    % Fitting an ohmic bath to a given set of decoherence times

    iTd = 1/T2 -0.5/T1; % inverse pure dephasing time
    if (iTd < 0)
      error('Unphysical decoherence times!')
    end
    
    % match bath couplings to T1, T2
    temp = b.scale * delta/2;

    alpha = atan2(1, sqrt(T1*iTd * coth(temp) * temp * b.cut_func(delta)));
    % dimensionless system-bath coupling factor squared
    N = iTd*b.scale/(4*pi*cos(alpha)^2);

    % qubit Hamiltonian
    H = -delta/2 * qit.sz;

    % noise coupling
    D = sqrt(N)*(cos(alpha)*qit.sz +sin(alpha)*qit.sx);

    % decoherence times in scaled time units
    %T1 = 1/(N * sin(alpha)^2 * 2*pi * delta * coth(temp) * b.cut_func(delta))
    %T_dephase = b.scale/(N *4*pi*cos(alpha)^2);
    %T2 = 1/(0.5/T1 +1/T_dephase)

  otherwise
    error('Unknown bath type.')
end
end
