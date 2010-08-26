function [p] = shor_factorization(N, cheat)
% SHOR_FACTORIZATION  Shor's factorization algorithm demo.
%  p = shor_factorization(N)       % simulate the full algorithm
%  p = shor_factorization(N, true) % cheat, avoid the quantum part
%
%  Simulates Shor's factorization algorithm, tries to factorize the integer N.
%
%  NOTE: This is a very computationally intensive quantum algorithm
%  to simulate classically, and probably will not run for any
%  nontrivial value of N (unless you choose to cheat, in which case
%  instead of simulating the quantum part we use a more efficient
%  classical algorithm for the order-finding).

%! P.W. Shor, "Algorithms For Quantum Computation: Discrete Logs and Factoring", Proc. 35th Symp. on the Foundations of Comp. Sci., 124 (1994).
%! M.A. Nielsen, I.L. Chuang, "Quantum Computation and Quantum Information" (2000), chapter 5.3.
% Ville Bergholm 2010


fprintf('\n\n=== Shor''s factorization algorithm ===\n\n')

global qit;

if (nargin < 2)
  cheat = false;
end
if (cheat)
  fprintf('(cheating)\n\n');
end

% number of bits needed to represent mod N arithmetic:
m = ceil(log2(N));
fprintf('Trying to factor N = %d (%d bits).\n', N, m);

% maximum allowed failure probability for the quantum order-finding part
epsilon = 0.25;

% number of index qubits required for the phase estimation
t = 2*m +1 +ceil(log2(2+1/(2*epsilon)));

fprintf('The quantum order-finding subroutine will need %d+%d qubits.\n\n', t, m);


% classical reduction of factoring to order-finding
while (true)
  a = 1+randi(N-2); % random integer, 1 < a < N
  fprintf('Random integer: a = %d\n', a)
  
  p = gcd(a, N);
  if (p ~= 1)
    % a and N have a nontrivial common factor p.
    % This becomes extremely unlikely as N grows.
    fprintf('Lucky guess, we found a common factor!\n')
    break
  end

  fprintf('Trying to find the period of f(x) = a^x mod N');

  if (cheat)
    % classical cheating shortcut
    r = find_order_cheat(a, N);
  else
    while (true)
      fprintf('.');
      % ==== quantum part of the algorithm ====
      [s1, r1] = find_order(a, N, t, m);
      [s2, r2] = find_order(a, N, t, m);
      % =============  ends here  =============

      if (gcd(s1, s2) == 1)
        % no common factors
        r = lcm(r1, r2);
        break;
      end
    end
  end
  fprintf('\n  =>  r = %d\n', r)

  % if r is odd, try again
  if (gcd(r, 2) == 2)
    % r is even

    x = mod_pow(a, r/2, N);
    if (mod(x, N) ~= N-1)
      % factor found
      p = gcd(x-1, N); % ok?
      if (p == 1)
        p = gcd(x+1, N); % no, try this
      end
      break
    else
      fprintf('a^(r/2) = -1 (mod N), try again...\n\n');
    end
  else
    fprintf('r is odd, try again...\n\n');
  end
end

fprintf('\nFactor found: %d\n', p)
end


function [s, r] = find_order(a, N, t, m)
% Quantum order-finding subroutine.
% Finds the period of the function f(x) = a^x mod N.

T = 2^t; % index register dimension
M = 2^m; % state register dimension

% applying f(x) is equivalent to the sequence of controlled modular multiplications in phase estimation
U = gate.mod_mul(a, N, M);

% state register initialized in the state |1>
st = state(1, M);

% run the phase estimation algorithm
reg = phase_estimation(t, U, st, true); % use implicit measurement to save memory

% measure index register
[dummy, res] = measure(reg, 1);

num = res-1;

% another classical part
r = find_denominator(num, T, T+1);
s = num*r/T;
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
  %if (temp/y < 1 / (2*max_den^2))
    return
  end

  % invert the remainder (swap numerator and denominator)
  x = y;
  y = temp;
end
end


function r = find_order_cheat(a, N)
% Classical order-finding algorithm.

for r = 1:N
  if (mod_pow(a, r, N) == 1)
    return
  end
end
end


function res = mod_pow(a, x, N)
% Computes mod(a^x, N) using repeated squaring mod N.
% x must be a positive integer.

b = fliplr(dec2bin(x) - '0'); % exponent in little-endian binary
res = 1;
for k = 1:length(b)
  if (b(k))
    res = mod(res*a, N);
  end
  a = mod(a*a, N); % square it
end
end
