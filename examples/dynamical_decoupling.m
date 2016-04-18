function dynamical_decoupling(seqs, titles)
% DYNAMICAL_DECOUPLING  Dynamical decoupling sequence demo.
%
%  Compare the performance of various dynamical decoupling
%  sequences against dephasing.
%

% Ville Bergholm 2016


fprintf('\n\n=== NMR dynamical decoupling sequences ===\n')

global qit
sx = qit.sx;
sy = qit.sy;
sz = qit.sz;
I = eye(2);


drive_strength = 4;

t = 10 * drive_strength;


if nargin < 1
    seqs = {seq.dd('wait', t), seq.dd('cpmg', t), seq.dd('uhrig', t, 3)};
    titles = {'Wait', 'Carr-Purcell-Meiboom-Gill', 'Uhrig-3'};
else
    if nargin < 2
        titles = {'User-given seq'};
    end
end





psi = state(kron(qit.p0, I/2), [2, 2]); % initial state

ns = length(seqs);
%fid = zeros(ns, ng, nf);
figure;

for q=1:ns
    subplot(1,ns,q)
    s = seqs{q};
    s.A = -1i*0.5*kron(sz, sz) * 0.1 / drive_strength;
    s.B = {-1i*0.5*kron(sx, I), -1i*0.5*kron(sy, I)};

    % add initial pi/2 rotation
    s.tau = [pi/2; s.tau];
    s.control = [0, 1; s.control];

    % apply sequence on state psi, plot the evolution
    [out, t] = seq.propagate(psi, s, @(s) s.ev(kron(sx, I)));

    plot(t, cell2mat(out));
    xlabel('t');
    ylabel('<y>');
    title(sprintf('%s, dephasing, rabi = %g', titles{q}, drive_strength));
    grid on;
end
