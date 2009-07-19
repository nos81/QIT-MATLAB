function out = adiabatic_propagate(s, H0, H1, t, steps)
% STATE/ADIABATIC_PROPAGATE  Propagate a quantum state in time adiabatically.
%  output = adiabatic_propagate(s, H0, H1, t, steps)
%
%  The Hamiltonian of the system is linearly interpolated between H0 and H1.

% Jacob D. Biamonte 2008
% Ville Bergholm 2009


if (nargin < 5)
  steps = 50;

  if (nargin < 4)
    error('Needs a state, H0, H1 and a time duration.');
  end
end

if (size(H0) ~= size(H1))
  error('Dimensions of the Hamiltonians do not match.');
end

d = size(H0, 2);
if (size(s.data, 1) ~= d)
  error('Dimension of the Hamiltonian does not match the dimension of the state.');
end

dt = t/steps;
dq = 1/steps;

[v,d] = eig(H1);
[S, I] = sort(diag(d), 'ascend');
lowest = v(:, I(1:4));

for k=1:steps
  out{k} = s; % store state

  q = (k-1)*dq;  % interpolation coefficient
  H = (1-q)*H0 + q*H1;
  U = expm(-i*H*dt);

  energies(:,k) = sort(eig(H), 'ascend'); % store energies

  if (size(s.data, 2) == 1)
    % state vector
    inner(:,k) = abs(s.data'*lowest).^2; % overlap with lowest final states

    s.data = U*s.data;    
  else
    % state operator
    inner(:,k) = diag(lowest'*s.data*lowest); % overlap with lowest final states

    s.data = U*s.data*U';
  end
end

out{k+1} = s;
energies(:,k+1) = sort(eig(H1), 'ascend'); % store energies
inner(:,k+1) = abs(s.data'*lowest).^2; % overlap with lowest final states

tt = 0:dt:t;


figure(1);
hold off;
plot(tt/t, energies);
%hold on;
grid on;
%m = 4;
%set(gca, 'XTick', [0:0.25*t:t]);
%set(gca, 'XTickLabel', {'s = 0';'s = 0.25';'s = 0.5';'s = 0.75'; 's = 1'});
%set(gca, 'YTick', [min(min(energies)):0.25*m:max(max(energies))])
%set(gca,'YTickLabel',{'0';'0.5';'1';'1.5'; '2'})
title('Energy spectrum');
xlabel('Adiabatic time');
ylabel('Energy');
axis([0, 1, min(min(energies)), max(max(energies))]);



figure(2);
hold off;
plot(tt/t, inner); %, 'LineWidth', 1.7);
%hold on;
grid on;
%set(gca, 'XTick', [0:0.25*t:t]);
%set(gca, 'XTickLabel', {'s = 0';'s = 0.25';'s = 0.5';'s = 0.75'; 's = 1'});
%set(gca, 'YTick', [0:0.5:1.1]);
%set(gca, 'YTickLabel', {'0';'0.5'; '1'});
title('Squared overlap of current state and H1 eigenstates');
xlabel('Adiabatic time');
ylabel('Probability');
%ylabel('|<l(s)|\psi(s)>|^2');
legend('\psi_0','\psi_1','\psi_2','\psi_3')
%axis([0, 1, 0, 1]);
axis([0, 1, 0, max(max(inner))]);



