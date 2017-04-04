function [res, s] = test_twist(a, sig)
% Produces different kinds of twist plots.

% Ville Bergholm 2011


global qit

if (false)
% random three-qubit correlation matrix
r = rand(4,4,4)-0.5; r(1,1,1) = 1;
r(2:4, 2:4, 2:4) = zeros(3,3,3); % no genuine 3q correlations
s = bloch_state(r); % into state op
s.data = s.data + eye(8); % make it positive
s = s/trace(s);
eig(s.data)
res = twist(s);

return
end

cluster = normalize(state('000')-state('001')-state('010')-state('011')+state('100')-state('101')+state('110')+state('111'));


x = state('001');
y = state('010');
z = state('100');
s0 = state('000');
s1 = state('111');

ghz = state('ghz');

sg = state('bell4'); % singlet
I = state(eye(2)/2, 2);
temp = tensor(sg, I);
b1 = normalize(temp +reorder(temp, [2 3]) +reorder(temp, [1 3]));


np = 30
p = linspace(0.05, 0.95, np);
nq = 30
q = linspace(0.0, 1.45, nq);
res = zeros(np, nq);

flip = mkron(qit.sx, qit.sx, qit.sx);

%sig = [1 1 1];

for j = 1:nq
  for k = 1:np
    w = normalize(a(1)*x +a(2)*y +a(3)*z +q(j)*s0);
    wbar = w.prop(flip);
    s = p(k)*to_op(w) +(1-p(k))*to_op(wbar);

    %SLOCC = mkron(rand_SU(2), rand_SU(2), rand_SU(2));
    %s = s.prop(SLOCC);
    
    %s = p(k)*to_op(ghz) +(1-p(k))*to_op(w);

    %s = p(k)*to_op(cluster) +(1-p(k))*to_op(b1);
    %s.data = (1-q(j))*s.data +q(j)*diag(diag(s.data));
    
    res(k,j) = twist(s, 1:3, sig);
  end
  j
end

norm(imag(res))
res = real(res);

figure
plot_pcolor(res, q, p, [0, max(max(res))]);
%pcolor(a, b, W);
%surf(q, p, res);
axis equal tight;
title(horzcat('twist, sig = [', sprintf('%d ', sig), '], w/wbar mixture, a = [', sprintf('%4f ', a), ']'))
%title(horzcat('twist, cluster/b1 mixture'))
xlabel('q')
ylabel('p')
colormap('default')
shading flat 

return


Sigma = [a 0 0 b; 0 d 0 0; 0 0 -d 0; c 0 0 a+c-b];

s = bloch_state(Sigma)
trace(s)
eig(s.data)
end
