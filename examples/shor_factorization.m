function [p] = shor_factorization(N)
% SHOR_FACTORIZATION  Shor factorization algorithm demo.
%  p = shor_factorization(N)
%
%  Simulate Shor factorization, try to factorize the integer N = p*q.
%
%  NOTE: This is a very computationally intense quantum algorithm to
%  simulate classically, and probably will not run for any
%  nontrivial value of N.

%!
% Ville Bergholm 2010


fprintf('\n\n=== Shor factorization algorithm ===\n')

global qit;

% number of bits needed to represent N
m = ceil(log2(N+1));
fprintf('Trying to factor %d (%d bits)\n\n', N, m);


% classical reduction to order-finding
while (true)
  a = 1+randi(N-2); % random integer, 1 < a < N
  fprintf('Random integer: a = %d\n', a)
  
  p = gcd(a, N);
  if (p ~= 1)
    % a and N have a nontrivial common factor p.
    % This becomes extremely unlikely as N grows.
    fprintf('Lucky guess!\n')
    break
  end

  r = find_order_cl(a, N, m); % quantum part of the algorithm
  fprintf('r = %d\n', r)

  % if r is odd, try again
  if (gcd(r, 2) == 2)
    % r is even
    temp = a^(r/2);
    if (mod(temp, N) ~= N-1)
      % factor found
      p = gcd(temp-1, N); % ok?
      if (p == 1)
        p = gcd(temp+1, N); % no, try this
      end
      break
    end    
  end
  fprintf('failed...\n\n')
end


fprintf('Factor found: %d\n', p)
end



function r = find_order_cl(a, N, m)
% classical algorithm
for r = 1:N
  if (mod(a^r, N) == 1)
    return
  end
end
end

function r = find_order(a, N, m)
% quantum order-finding subroutine
% finds the period of the function f(x) = a^x mod N

epsilon = 0.2; % failure probability

% number of index qubits required for the phase estimation
t = 2*m +1 +ceil(log2(2+1/(2*epsilon)));

fprintf('Quantum order-finding subroutine: Using %d+%d qubits\n', t, m);
fprintf('Trying to find the period of f(x) = %d^x mod %d  =>  ', a, N);


%dim = [2*ones(1, t), 2*ones(1, m)];
T = 2^t;
L = 2^m;

% index register initialized in uniform superposition (skip the Hadamards)
r1 = state(ones(T, 1)/sqrt(T), T);
% Hadamard gates, initial superposition
%H = gate.walsh(t);
%r1 = u_propagate(r1, H);

% state register initialized in state |1>
r2 = state(1, L);

reg = tensor(r1, r2);

% apply f(x) (equivalent to the sequence of controlled unitaries in phase estimation)
U = gate.mod_mul(a, N, L);
for k = 1:t
  ctrl = -ones(1, t);
  ctrl(k) = 1;
  temp = gate.controlled(U^(2^(t-k)), ctrl);
    
  reg = u_propagate(reg, temp);
end

% do an inverse quantum Fourier transform
QFT = gate.qft(t);

%reg_t = ptrace(reg, [2]); % trace over state register
%reg = u_propagate(reg_t, QFT');

reg = u_propagate(reg, kron(QFT', speye(L)));


% measure first register
[~, res] = measure(reg, 1);

% another classical part
r = find_denominator(res-1, T, T+1);
end


function a = continued_fraction(p, q)
% Returns the vector of quotients a for the continued fraction representing x = p/q.

a = [];
while (1)
  a(end+1) = floor(p/q); % integer part

  temp = mod(p, q); % remainder
  if (temp == 0)
    break;
  end

  % invert
  p = q;
  q = temp;
end
end

function d_1 = find_denominator(x, y, max_den)
% Finds the denominator q for p/q \approx x/y such that q < max_den
% using a continued fraction representation for x.
%
% We use floor and mod here, which could be both implemented using
% the classical Euclidean algorithm.

d_2 = 1;
d_1 = 0;

while (1)
  a = floor(x/y); % integer part == a_n
  temp = a*d_1 +d_2; % n:th convergent denumerator d_n = a_n*d_{n-1} +d_{n-2}

  if (temp >= max_den)
    return % stop
  end
  d_2 = d_1;
  d_1 = temp;
 
  temp = mod(x,y); % x - a*y; % subtract integer part

  if (temp == 0)
  %if (temp/y < 0.5 / max_den^2)
    return
  end

  % invert the remainder (swap numerator and denominator)
  x = y;
  y = temp;
end
end
