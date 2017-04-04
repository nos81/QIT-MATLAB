function s = prop(s, P)
% PROP  Propagate the state a finite step in time.
%  q = prop(s, P)
%
%  Propagates the state s using the propagator P, which is either a
%  unitary Hilbert space propagator, a general Liouville space propagator,
%  or a cell vector of Kraus operators.
%
%  Returns the resulting state.

% TODO allow the user to apply P only to some subsystems of s (pad with identity)

% Ville Bergholm 2009-2017

if nargin < 2
    error('Propagator required.');
end

if ~iscell(P)
    P = {P};
end
n = length(P);

%temp = 0;
%for k=1:n
%  temp = temp + P{k}'*P{k};
%end
%if (norm(temp - eye(size(temp))) > qit.tol)
%  warning('Unphysical quantum operation.')
%end
% TODO: If n > prod(s.dims())^2, there is a simpler equivalent
% quantum operation. Should the user be notified?

% are they lmaps?
if isa(P{1}, 'lmap')
    % is the input dimension correct?
    if ~is_concatenable(P{1}, s)
        error('zzz')
    end
    dim = P{1}.dim{1}; % output dimension
    if s.is_ket()
        s.dim = {dim, 1};
    else
        s.dim = {dim, dim};
    end
    for k=1:n
        P{k} = P{k}.data;
    end
end

if n == 1
    % a unitary or a Liouville propagator
    if length(P{1}) == size(s.data, 1)
        % unitary
        if is_ket(s)
            % state vector
            s.data = P{1} * s.data;
        else
            % state operator
            s.data = P{1} * s.data * P{1}';
        end
    else
        % Liouvillian propagator
        s = to_op(s);
        s.data = inv_vec(P{1} * vec(s.data));
    end
else
    % set of Kraus operators
    s = to_op(s); % into a state operator
    temp = 0;
    for k=1:n
        temp = temp +P{k} * s.data * P{k}';
    end
    s.data = temp;
end
