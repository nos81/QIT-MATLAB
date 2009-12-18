function [out, t] = seq_propagate(s, seq, out_func)
% STATE/PROPAGATE  Propagate the state in time using a control sequence.
%  [out, t] = propagate(s, seq, out_func)
    
% Ville Bergholm 2009


global qit;

if (nargin < 3)
    out_func = @(x) x; % no output function given, use a NOP
    
    if (nargin < 2)
        error('Needs a stuff');
    end
end

if (s.dim(1) ~= 2)
  error('only works on qubits for now.');
end

base_dt = 0.1;
n = size(seq, 1); % number of pulses
t = [0];
out{1} = out_func(s);

for k=1:n
    H = 0.5*(qit.sx*seq(k, 1) +qit.sy*seq(k, 2) +qit.sz*seq(k, 3));
    T = seq(k, end);
    
    n_steps = ceil(T/base_dt);
    dt = T/n_steps;

    U = expm(-i*H*dt);
    for j=1:n_steps
        s = u_propagate(s, U);
        out{end+1} = out_func(s);
    end

    temp = t(end);
    t = [t, linspace(temp+dt, temp+T, n_steps)];
end
