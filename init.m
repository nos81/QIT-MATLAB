% Initialization script for the quantum information toolkit.

addpath(pwd);
addpath([pwd '/utils']);
addpath([pwd '/examples']);

% define global variables in the toolkit namespace
global qit;

% toolkit version number
qit.version = '0.9.10';
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

% some relevant physical constants (CODATA 2010)
qit.hbar = 1.054571726e-34; % Planck constant / (2 pi), J s
qit.kB   = 1.3806488e-23;   % Boltzmann constant, J/K
qit.eV   = 1.602176565e-19; % electron volt, J
qit.c    = 299792458;       % speed of light in vacuum, m/s


% caching
qit.gellmann = {};
qit.tensorbasis = {};
qit.tensorbasislocal = {};
qit.angular_momentum = {};
qit.qft = {};
