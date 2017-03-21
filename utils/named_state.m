function rho = named_state(name, dim, p)
% NAMED_STATE  Construct different types of named quantum states.
%  rho = named_state(name, dim, params)
%
%  The currently supported named states are
%    bell1, bell2, bell3, bell4
%    ghz         (Greenberger-Horne-Zeilinger states)
%    w           (W states)
%    werner
%    isotropic

% Ville Bergholm 2017


name = lower(name); % convert to lower case

if nargin < 2
    dim = [2, 2];
end
switch name
  case {'ghz', 'w'}
    if nargin == 1
        dim = [2, 2, 2];
    end
end

n = length(dim);
s = zeros(prod(dim), 1);

switch name
  case {'bell1', 'bell2', 'bell3', 'bell4'}
    % 2-qubit Bell state
    if ~isequal(dim, [2, 2])
        error('Standard Bell states must consist of two qubits.')
    end
    Q_Bell = [0 0 1 1i; 1 1i 0 0; -1 1i 0 0; 0 0 1 -1i] / sqrt(2);
    rho = state(Q_Bell(:, name(5)-'0'), dim);

  case 'ghz'
    % Greenberger-Horne-Zeilinger state
    s(1) = 1;
    s(end) = 1;
    rho = state(s, dim).normalize();

  case 'w'
    % W states
    ind = 1;
    for k=n:-1:1
        s(ind*(dim(k) - 1) + 1) = 1; % MATLAB indexing starts at 1
        ind = ind*dim(k);
    end
    rho = state(s, dim).normalize();

  case 'werner'
    rho = state(gate.werner(p, dim));
  case 'isotropic'
    rho = state(gate.isotropic(p, dim));
  otherwise
    error('Unknown named state ''%s''.', name)
end
