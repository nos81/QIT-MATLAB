function [out, t] = seq_propagate(s, seq, out_func)
% SEQ_PROPAGATE  Propagate the state in time using a control sequence.
%  [out, t] = propagate(s, seq, out_func)
    
% Ville Bergholm 2009-2011


if (nargin < 3)
    out_func = @(x) x; % no output function given, use a NOP
    
    if (nargin < 2)
        error('Needs a control sequence.');
    end
end


base_dt = 0.1;
n = length(seq.tau); % number of pulses
t = [0];
out{1} = out_func(s);

for j=1:n
    G = seq.A;
    for k = 1:length(seq.B)
        u = seq.control(j, k);
        G = G + u * seq.B{k};
    end

    T = seq.tau(j);
    n_steps = ceil(T / base_dt);
    dt = T / n_steps;

    P = expm(-dt * G);
    for q=1:n_steps
        s = u_propagate(s, P);
        out{end+1} = out_func(s);
    end

    temp = t(end);
    t = [t, linspace(temp+dt, temp+T, n_steps)];
end
