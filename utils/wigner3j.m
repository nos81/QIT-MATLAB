function ret = wigner3j(a, b, c, x, y, z)
% WIGNER3J  Wigner 3-j symbols.
%  wigner3j(a, b, c, x, y, z)
%
%  Returns the Wigner 3-j coefficients used in angular momentum coupling.
%  Uses the Racah formula.

% Ville Bergholm 2017


% triangle condition for the spins
if abs(a-b) > c || a+b < c
    ret = 0;
    return
end
% z component condition
if abs(x) > a || abs(y) > b || abs(z) > c
    ret = 0;
    return
end
% triangle condition for the z components
if x+y+z ~= 0
    ret = 0;
    return
end

% triangle coefficient
ret = factorial(-a+b+c)*factorial(a-b+c)*factorial(a+b-c)/factorial(a+b+c+1);
ret = ret * factorial(a+x)*factorial(a-x)*factorial(b+y)*factorial(b-y)*factorial(c+z)*factorial(c-z);

temp = 0;
tmin = max([0, b-c-x, a-c+y])
tmax = min([-c+a+b, a-x, b+y])
for t=tmin:tmax
    t
    temp = temp +(-1)^t / (factorial(t)*factorial(c-b+x+t)*factorial(c-a-y+t)*factorial(a+b-c-t)*factorial(a-x-t)*factorial(b+y-t));
end
ret = (-1)^(a-b-z) * sqrt(ret) * temp;
