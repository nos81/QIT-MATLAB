function s = cpmg(t, n)
% CPMG  Carr-Purcell-Meiboom-Gill sequence.
%  s = cpmg(t, n)
%
%  Returns the Carr-Purcell-Meiboom-Gill sequence of n repeats with waiting time t.
%  The purpose of the CPMG sequence is to facilitate a T_2 measurement
%  under a nonuniform z drift, it is not meant to be a full memory protocol.
%  The target operation for this sequence is identity.
    
% Ville Bergholm 2007-2011


s = seq.nmr([pi/2, pi/2]); % initial y rotation

% step: wait, pi x rotation, wait
step_tau  = [t; pi; t];
step_ctrl = [0 0; 1 0; 0 0]; 

for k=1:n
  s.tau = [s.tau; step_tau];
  s.control = [s.control; step_ctrl];
end
