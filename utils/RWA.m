function [Q, Q_slow, omega_rot, a_min, a_max] = RWA(H0, A, carrier, verbosity, allow_crosstalk, safety_factor)
% RWA  Rotating wave approximation.
%  [Q_static, Q_slow, omega_rot, a_min, a_max] = RWA(H0, A [, omega_carrier])
%
%  Performs the rotating frame transformation U = \exp(-i H0 t), where
%  H0 is the hermitian generator of the transformation, on the
%  hermitian operator A, giving A'(t) = U^\dagger A U.
%  The fast-rotating terms are then thrown away.
%
%  The matrix Q_static is a sum of the (nearly) static terms of A'.
%
%  Q_slow is a cell vector of function handles @(t) representing the
%  terms rotating with a slow but nonvanishing rate.
%
%  The corresponding angular frequencies are given in the omega_rot vector.
%
%  a_min is a vector, each element denoting the
%  control amplitude above which the corresponding term in Q_slow
%  first passes the safety_factor limit and thus becomes significant.
%  The higher the a_min, the less significant the term.
%
%  If omega_carrier is given, A is assumed to oscillate in time,
%
%    A(t) = A a(t) \cos(\omega_carrier t +\phi(t)).
%
%  The maximum amplitude a_max = max_t |a(t)| is usually chosen low enough so
%  that only the transition(s) closest to omega_carrier are excited (no crosstalk),
%  and the fast-rotating terms can be neglected.
%  In this case A' = a(t) Q(t, phi).
%  Q_static is a function handle @(phi) and Q_slow a cell vector of handles @(t, phi).
%
%  The inputs are dimensionless (Hamiltonians divided by \hbar and multiplied by TU,
%  angular frequencies multiplied by TU).

% Ville Bergholm 2013-2017


% The logic here is the following:
%
% For each H0 eigenspace there is a corresponding eigenfrequency omega_k and a hermitian projector P_a.
% Each (a,b)-projected "block" A_ab = P_a*A*P_b rotates with the frequency delta_ab = omega_a -omega_b.
% cos(omega_carrier t) consists of two counter-rotating terms.
% Combined with delta_ab this yields a faster and a slower
% component, |omega_fast| >= |omega_slow| (equal if omega_carrier = 0 or a = b).
% The term exp(i omega_xxxx t) a_max A_ab can be neglected if the total rotation speed
% |omega_xxxx| > safety_factor * |A_ab| * a_max.
% To obtain a fixed A', we must be able to ignore all the faster components,
% and the slower components of unwanted transitions (crosstalk), which is
% done by setting a_max low enough.


% omega-tolerance for combining rotating terms/blocks/modes
tol_dH = 1e-3;
% |A|-tolerance for keeping block groups
tol_norm = 0;
% relative omega-tolerance for deciding which modes are unwanted
% and, if crosstalk is to be avoided, included in a_max determination
tol_target = 1e-7;
% ignored terms with |omega|/|A| ratios above this will not be even shown
tol_show = 1e4;
% |omega|-tolerance for considering terms static
tol_veryslow = 1e-3;

% formatting string for output
format = ' %s:  2pi*|%6.4g| / %6.4g = %6.4g\n';


if nargin < 6
    safety_factor = 500; % TODO justify

if nargin < 5
    allow_crosstalk = false;

if nargin < 4
    verbosity = 1;

% If carrier > 0 is given, H is multiplied by cos(carrier*t).
if nargin < 3
    carrier = 0;
elseif carrier <= 0
    error('Carrier frequency must be positive.')
end
end
end
end

% The identity component of a Hamiltonian A plays no physical role,
% but affects the op_norm calculation later on when a == b.
% Removing the id. comp. always lowers the Frobenius norm.
dim = length(A);
id_comp = trace(A) * eye(dim) / dim;
A = A -id_comp;

[dH, AA, group_label] = rotating_frame(H0, A,...
                                       'tol_dH',  tol_dH, ...
                                       'tol_norm',tol_norm);

%% compute volatilities

% compute norms
A_norm = zeros(size(dH));
for k=1:length(dH)
    % norm of the summed A projection
    A_norm(k) = op_norm(AA{k});
    % TODO with carrier > 0, drop diagonal blocks since they will
    % be discarded later anyway?
    if dH(k) ~= 0  % offdiagonal
        % the +h.c. terms counts towards the norm, too
        A_norm = A_norm * sqrt(2);
        if carrier ~= 0
            % slow or fast mode alone (for carrier == 0 they're the same)
            A_norm = A_norm * 0.5;
        end
    end
end
omega_fast = carrier +dH;
omega_slow = carrier -dH;
% All omega_slow <= carrier <= omega_fast, with equality only for diagonal blocks (dH = 0).
abs_omega_slow = abs(omega_slow);
% "volatility" of the term. The higher, the less important.
vol_slow = abs_omega_slow ./ A_norm;
vol_fast = abs(omega_fast) ./ A_norm;

if verbosity >= 2
    figure();
    semilogy(omega_fast/2/pi, A_norm, 'r+', abs_omega_slow/2/pi, A_norm, 'bo')
    hold on
    text(abs_omega_slow/2/pi, A_norm, group_label);
end



%% in the carrier mode, determine a_max such that the fast mode(s) (and maybe crosstalk) do not interfere

if carrier ~= 0
    % eliminate fast modes
    % NOTE On the diagonal |omega_slow| == |omega_fast|, so this
    % also eliminates diagonal blocks entirely when there is a carrier.
    % We thus cannot target the diagonal blocks... which would be crazy anyway?
    [min_vol_fast, j] = min(vol_fast);

    % find lowest |omega_slow|
    temp = min(abs_omega_slow);

    % all the slower modes within tol_target*carrier of the lowest one are the 'targeted' ones
    targets = find(abs_omega_slow <= temp +tol_target * carrier);
    unwanted = setdiff(1:length(abs_omega_slow), targets);

    ignore_crosstalk = allow_crosstalk || isempty(unwanted);
    if ignore_crosstalk
        min_vol_slow = inf; % do not care about crosstalk
    else
        % eliminate crosstalk: all slow modes except the targeted
        % ones must have negligible effect.
        [min_vol_slow, k] = min(vol_slow(unwanted));
    end

    if verbosity >= 1
        fprintf('\nf_carrier: %g, targets: ', carrier/2/pi);
        fprintf('%s, ', group_label{targets});
        fprintf('\n');
        if ~ignore_crosstalk
            fprintf('crosstalk limit %s: vol < %g\n', group_label{unwanted(k)}, min_vol_slow);
        end
	% safety factor FIXME not printed!!
        fprintf('fast mode limit %s: vol < %g\n', group_label{j}, min_vol_fast);
        %fprintf(['fast mode limit', format], group_label{j}, omega_fast(j)/2/pi, A_norm(j), min_vol_fast);
    end
    min_vol_discard = min(min_vol_fast, min_vol_slow);
    Q = @(phi) 0;
else
    % no carrier, no amplitude, no "unwanted" modes, discard decision is based on an absolute limit.
    min_vol_discard = safety_factor;
    if verbosity >= 1
        fprintf('\nf_carrier: %g\n', carrier/2/pi);
    end
    Q = 0;
end

if verbosity >= 2
    grid on
    a = axis();
    line([1, 1] * carrier/2/pi, [a(3), a(4)]);
    temp = linspace(1, a(2), 500);
    line(temp, temp * 2*pi/min_vol_discard);
    xlabel('f (1/TU)')
    ylabel('|P_a A P_b|')
    title(sprintf('f_{carrier}: %g, min_vol_discard: %g', carrier/2/pi, min_vol_discard))
end


%% see which slowly rotating terms to keep

Q_slow = {};
a_min = [];
omega_rot = [];

noshow = '';

% loop over the groups, in order of ascending vol_slow
[~, ind] = sort(vol_slow);
for k=ind.'
    if verbosity >= 1
        % print the label describing this group
        label = sprintf(format, group_label{k}, omega_slow(k)/(2*pi), A_norm(k), vol_slow(k));
    end

    if vol_slow(k) >= 0.9999 * min_vol_discard; % numerical safety
        % Even the slower mode is too fast, ignore.
        if verbosity >= 1
            if vol_slow(k) <= tol_show
                fprintf(['ignore  ', label])
            else
                noshow = [noshow, group_label{k}];
            end
        end
        continue
    end
    % omega_slow term is slow enough to keep.

    % build the rotating term
    if carrier == 0
        % slow and fast modes are the same, no phi
        R = @(t) exp(1i*omega_slow(k)*t) * AA{k};
    else
        % ignore the fast mode (min_vol_discard takes care of them), hence the 0.5.
        % diagonal modes have also been discarded.
        R = @(t, phi) 0.5 * exp(1i*(omega_slow(k)*t +phi)) * AA{k};
    end

    % rotating_frame only returns the "upper triangle" of projectors, so
    if dH(k) ~= 0
        % off-diagonal group, add the hermitian conjugate
        S = @(x) x+x';
    else
        % do nothing
        S = @(x) x;
    end

    if abs_omega_slow(k) <= tol_veryslow
        % Very slow component wrt. TU, add directly to Q.
        if verbosity >= 1
            fprintf(['fixed   ', label])
        end
        if carrier == 0
            Q = Q +S(R(0));
        else
            % HACK recursive function definition... evaluating
            % it and saving the results before use improves performance.
            Q = @(phi) Q(phi) +S(R(0, phi));
        end
    else
        % Slowly rotating component.
        % Can happen anywhere with carrier > 0, or with carrier = 0 in off-diagonal groups.
        % amplitude under which the volatility of the term exceeds the safety_factor
        a_min(end+1) = vol_slow(k) / safety_factor;
        % rotation speed
        omega_rot(end+1) = omega_slow(k);
        if verbosity >= 1
            fprintf(['rotating', label])
        end
        if carrier == 0
            Q_slow{end+1} = @(t) S(R(t));
        else
            Q_slow{end+1} = @(t, phi) S(R(t, phi));
        end
    end
end

% by driving with this amplitude, the first discarded mode has a
% volatility equal to safety_factor
a_max = min_vol_discard / safety_factor;

% TEST restore the identity component to a fixed term to make it look neater
if carrier == 0
    Q = Q +id_comp;
end

if verbosity >= 1
    if ~isempty(noshow)
        fprintf('ignore+  %s\n', noshow);
    end
    fprintf('safety: %g\n', safety_factor);
    fprintf('a_max: %g\n', a_max);
    fprintf('a_min:');
    fprintf(' %g,', a_min);
    fprintf('\nomega_rot:');
    fprintf(' %g,', omega_rot);
    fprintf('\n');
end
end


function ret = op_norm(A)
% Given a traceless Hermitian operator A, returns a positive scalar describing the
% magnitude of the unitary flow generated by it.
  ret = norm(A, 'fro');
end
