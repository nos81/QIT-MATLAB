function s = dd(name, t, n)
% DD  Dynamical decoupling and refocusing sequences.
%  s = dd(name, t, n, m)
%
%  name    name of the sequence: hahn, cpmg, uhrig, xy4
%  t       total waiting time
%  n       order (if applicable)

%  The target operation for these sequences is identity.

%! G.Uhrig, "Keeping a Quantum Bit Alive by Optimized \pi-Pulse Sequences", PRL 98, 100504 (2007).
% Ville Bergholm 2007-2016

% Multiplying the pi pulse strength by factor s is equivalent to A -> A/s, t -> t*s.


% which sequence?
switch name
  case 'wait'
    % Do nothing, just wait.
    tau = 1;
    phase = [];

  case 'hahn'
    % Basic Hahn spin echo
    tau = [0.5, 0.5];
    phase = [0];

  case 'cpmg'
    % Carr-Purcell-Meiboom-Gill
    % The purpose of the CPMG sequence is to facilitate a T_2 measurement
    % under a nonuniform z drift, it is not meant to be a full memory protocol.
    tau = [0.25, 0.5, 0.25];
    phase = [0, 0];

  case 'uhrig'
    % Uhrig's family of sequences
    % n=1: Hahn echo
    % n=2: CPMG
    delta = 0:n+1;
    delta = sin(pi * delta/(2*(n+1))).^2;
    tau = delta(2:n+2)-delta(1:n+1);  % wait durations
    phase = zeros(1, n);

  case 'xy4'
    % uncentered version
    tau = [1, 1, 1, 1] / 4;
    phase = [0, pi/2, 0, pi/2];
  otherwise
    error('Unknown sequence.')
end

% initialize the sequence struct
s = seq.nmr([0, 0]);

% waits and pi pulses with given phases
ind = 1;
for k=1:length(phase)
    % wait
    s.tau(ind,1) = t * tau(k);
    s.control(ind, :) = [0, 0];
    ind = ind+1;
    % pi pulse
    p = phase(k);
    s.tau(ind,1) = pi;
    s.control(ind, :) = [cos(p), sin(p)];
    ind = ind+1;
end

if length(tau) > length(phase)
    % final wait
    s.tau(ind,1) = t * tau(end);
    s.control(ind, :) = [0, 0];
end
