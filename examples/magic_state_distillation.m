function [rho] = magic_state_distillation(type)
% Magic state distillation demo.
%  E = magic_state_distillation(type)
%
%  dgfd

%@article{,
%  title = {Universal quantum computation with ideal Clifford gates and noisy ancillas},
%  author = {Bravyi, Sergey and Kitaev, Alexei},
%  journal = {Phys. Rev. A},
%  volume = {71},
%  issue = {2},
%  pages = {022316},
%  numpages = {14},
%  year = {2005},
%  doi = {10.1103/PhysRevA.71.022316}, }
% Ville Bergholm 2014


fprintf('\n\n=== Magic state distillation ===\n')

global qit;

% Clifford generators
H = qit.H;
K = diag([1, 1i]);

% Pauli matrices
I = eye(2);
X = qit.sx;
Y = qit.sy;
Z = qit.sz;


% random initial state
rho = state(rand_positive(2)).bloch_vector();

% into the positive octant (e.g. using Cliffords)
rho = bloch_state(abs(rho))


% what kind of magic state?
switch type
  case 'T'
    x = 1/sqrt(3);
    T0 = to_ket(bloch_state([1, x, x, x], 2));
    T1 = to_ket(bloch_state([1, -x, -x, -x], 2));

    % initial state error
    ep = 1 -rho.fidelity(T0)^2
    ep0 = 0.5 * (1 -sqrt(3/7))
    if ep >= ep0
        disp('Distillation will not work with this initial state.')
    end
    
    T = exp(1i*pi/4) * K * H;
    
    % dephase rho (does not affect the T0 fidelity)
    temp = (rho +rho.u_propagate(T) +rho.u_propagate(T')) / 3;
    
    % use five copies
    rho_in = temp.tensorpow(5);

    % 5-qubit stabilizer code
    S{1} = mkron(X, Z, Z, X, I);
    S{2} = mkron(I, X, Z, Z, X);
    S{3} = mkron(X, I, X, Z, Z);
    S{4} = mkron(Z, X, I, X, Z);
    
    d = 2*ones(1, 5);
    U = lmap(mkron(Z*H, H, H, Z*H, I), {d, d});
    U = gate.controlled(XZYI, [1]).reorder() * U;

    % orthogonal projector to codespace
    temp = eye(2^5);
    PP = temp / 16;
    for k=1:4
        PP = (temp +S{k}) * PP;
    end

    % trivial syndrome state
    rho_s = rho_in.u_propagate(PP);
    % and the probability of obtaining it
    p_s = real(rho_s.trace())
    
    % decoding
    rho_out = rho_s.u_propagate(U_dec);
    
  case 'H'
    target = bloch_state([1, 0, 1] / sqrt(2), 2)  %.to_ket()
    
  otherwise
    error('Unknown state.')
end


