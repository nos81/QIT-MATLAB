function p = bernstein_vazirani(n, linear)
% BERNSTEIN_VAZIRANI  Bernstein-Vazirani algorithm demo.
%  p = bernstein_vazirani(n, linear)
%
%  Simulates the Bernstein-Vazirani algorithm, which, given a black box oracle
%  implementing a linear Boolean function f_a(x) := a \cdot x, returns the bit
%  vector a (and thus identifies the function) with just a single oracle call.
%  If the oracle function is not linear, the algorithm will fail.

%! See \cite{BV}.
% Ville Bergholm 2011-2014

fprintf('\n\n=== Bernstein-Vazirani algorithm ===\n')
fprintf('\nUsing %d qubits.', n);

if nargin < 2
    linear = true;
end


dim = qubits(n);
H = gate.walsh(n); % n-qubit Walsh-Hadamard gate
N = 2^n;

% random black box oracle encoding the Boolean function f (given as the diagonal)
if linear
    % linear f
    a = (rand(1, n) > 0.5) +0;
    f = linear_func(a);
    fprintf('\nUsing the linear function f_a(x) := dot(a, x), defined by the binary vector\n');
    a
else
    % general f
    f = (rand(1, N) > 0.5) +0;
    % special case: not(linear)
    %a = (rand(n, 1) > 0.5) +0; 
    %f = 1-linear_func(a);
    fprintf('\nNonlinear function f:\n');
    f
end
U_oracle = oracle(f);

% start with all-zero state
s = state(0, dim);
% initial superposition
s = s.u_propagate(H);
% oracle phase flip
s.data = U_oracle .* s.data;
% final Hadamards
s = s.u_propagate(H);

figure();
s.plot();
tt = 'Bernstein-Vazirani algorithm, ';
if linear
    temp = 'linear oracle';
else
    temp = 'nonlinear oracle (fails)';
end
tt = [tt, temp];
title(tt);

[p, res] = measure(s);
% measured binary vector
b = unravel_index(res, dim) -1;
fprintf('\nMeasured binary vector\n');
b
if ~linear
    g = linear_func(b);
    fprintf('\nCorresponding linear function g_b:\n');
    g
    fprintf('\nNormalized Hamming distance: |f - g_b| = %g.\n', sum(abs(f-g)) / N);
end
end


function U = oracle(f)
% Returns a unitary oracle for the Boolean function f(x), given as a truth table.
    U = (-1) .^ f.';
end
    
function f = linear_func(a)
% Builds the linear Boolean function f(x) = a \cdot x as a truth table.
    n = length(a);
    dim = qubits(n);
    N = prod(dim);
    f = zeros(1, N);
    for ind = 1:N
        % convert each index ind into the corresponding bitstring x
        x = unravel_index(ind, dim) -1;  % into zero-based multi-index
        % compute f(x)
        f(ind) = mod(dot(a, x), 2);
    end
end