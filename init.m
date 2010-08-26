% Initialization script for the quantum information toolkit.

addpath(pwd);
addpath([pwd '/utils']);
addpath([pwd '/examples']);

% define global variables in the toolkit namespace
global qit;

% toolkit version number
qit.version = '0.9.8';
fprintf('MATLAB Quantum Information Toolkit, version %s\n', qit.version)

% Pauli matrices
qit.I  = eye(2);
qit.sx = [0 1; 1 0];
qit.sy = [0 -i; i 0];
qit.sz = [1 0; 0 -1];

% qubit projectors
qit.p0 = [1 0; 0 0];
qit.p1 = [0 0; 0 1];

% easy Hadamard
qit.H = [1 1; 1 -1]/sqrt(2);

% error tolerance
qit.tol = max(1e-10, eps);

% some relevant physical constants (CODATA 2006)
qit.hbar = 1.054571628e-34; % J s
qit.kB   = 1.3806504e-23;   % J/K
qit.eV   = 1.602176487e-19; % J

% caching
qit.gellmann = {};
qit.tensorbasis = {};
qit.angular_momentum = {};
qit.qft = {};
