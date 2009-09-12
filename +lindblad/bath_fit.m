function [b, H, D] = bath_fit(b, delta, T1, T2)
% LINDBLAD/BATH_FIT  Setup a bath that reproduces given decoherence times for a qubit.
%
%   [bath, H, D] = bath_fit(bath, delta, T1, T2)
%
%  Sets the bath parameters to values which reproduce the required
%  decoherence times T1 and T2 (in s) for a single-qubit system.
%  delta is the energy splitting for the qubit (in units of hbar*omega0).
%
%  The function also returns the corresponding system Hamiltonian H and
%  the system-bath coupling operator D.

% Ville Bergholm 2009


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
    b.N = iTd*b.scale/(b.omega0*4*pi*cos(alpha)^2);

    % qubit Hamiltonian
    H = -delta/2 * qit.sz;

    % noise coupling
    D = cos(alpha)*qit.sz +sin(alpha)*qit.sx;

    % decoherence times in scaled time units
    T1 = 1/(b.N * sin(alpha)^2 * 2*pi * delta * coth(temp) * b.cut_func(delta))
    T_dephase = b.scale/(b.N *4*pi*cos(alpha)^2);
    T2 = 1/(0.5/T1 +1/T_dephase)

  otherwise
    error('Unknown bath type.')
end
end
