function s = cpmg(t, n)
% CPMG  Carr-Purcell-Meiboom-Gill sequence of n repeats with waiting time t.
%  s = cpmg(t, n)
%
%  The purpose of the CPMG sequence is to facilitate a T_2 measurement
%  under a nonuniform z drift, it is not meant to be a full memory protocol.
%  The target operation for this sequence is identity.
    
% Ville Bergholm 2007-2009


s = seq.nmr([pi/2, pi/2]); % initial y rotation

step = [0 0 0 t; seq.nmr([pi, 0]); 0 0 0 t]; % wait, pi x rotation, wait

for k=1:n
  s = [s; step];
end
