function [V, Sigma, W, S, errflag] = lorentz_svd(s)
% LORENTZ_SVD  Lorentz decomposition of a two-qubit state.
%
%  [V, Sigma, W, S] = lorentz_svd(s)
%
%  Decomposes the correlation matrix S of the two-qubit state s as
%  S = V * \Sigma * W.', where V and W are in SO^+(1, 3) and
%  \Sigma is diagonal.
%
%  Cannot yet handle all the possible cases.
    
%! F. Verstraete et al., "Local filtering operations on two qubits", PRA 64, 010101 (2001). 
% Ville Bergholm 2011


global qit eta T tol errflag

errflag = false;
tol = qit.tol * 1e3; % jordan isn't too precise

eta = diag([1, -1, -1, -1]);
T   = diag([-1, 1, 1, 1]);

% test: "infinitesimal" fudge factor to avoid numerical problems
% no good, twist isn't continuous!
%epsilon = 0;
%s = (1-epsilon)*to_op(s) +epsilon * state(rand_positive(4), [2 2]);

S = bloch_vector(s);
% eliminate possible imaginary errors
if (norm(imag(S)) > tol)
  error('imaginary stuff in S');
end
S = real(S);

% construct V
temp = S*eta*S.'*eta;
[V, Sigma] = xxx(temp, S);

% now W
temp = S.'*eta*S*eta;
[W, Sigmaw] = xxx(temp, S);

if (norm(Sigma-Sigmaw) > tol)
  error('sigmas don''t match')
end

% restore Sigma's signature
if (det(S) < 0)
  Sigma = Sigma * eta; % other possibility
end

A = (V\S)/W.';
% patch up holes in A with anything invertible (identity)
warn = false;
for k=1:4
    if (abs(Sigma(k,k)) < tol)
       A(k,k) = 1;
       warn = true;
    else
       A(:,k) = A(:,k)/Sigma(k,k);
    end
end

if (warn)
    disp('Warning: Sigma not invertible, parallel transporter ill-defined.')
end

V = V*A;

%norm(V*Sigma*W.' -S)
return

% SO^+ sign subgroup, commutes with diagonal Sigma, so we only need to fix V
Q{1} = eye(4);
Q{2} = diag([1 1 -1 -1]);
Q{3} = diag([1 -1 -1 1]);
Q{4} = diag([1 -1 1 -1]);



for k=1:4
  if (norm(V*Q{k}*Sigma*W.' -S) < tol)
    % found a working signature
    V = V * Q{k};
    return
  end
end

disp('no working signature found.')
V = NaN;
end


function [V, Sigma] = xxx(C, S)

global tol eta T errflag

if (rank(C) < rank(S))
  disp('deficient rank case')
  errflag = true;
  % TODO limited case (3)
end

% S = V Sigma W^T
% S eta S^T eta = V (Sigma eta Sigma^T eta) V^{-1}
% S^T eta S eta = W (Sigma^T eta Sigma eta) W^{-1}

[V, Sigma2] = eigsort(C);
% jordan would be more appropriate, the matrix is not always diagonalizable
% however that's a zero measure subset so...
% here we implicitly assume case (1) from the ref. (Sigma is diagonal)
% TODO case (3)

if (norm(V*Sigma2/V -C) > tol)
  disp('eig failed, need jordan')
  errflag = true;
end

% A takes care of the freedom in the eigendecomposition
% [Sigma^2, A] = 0, A invertible  =>
% C = V Sigma^2 V^{-1} = V A^{-1} Sigma^2 (V A^{-1})^{-1}

% now make sure V A^{-1} is an SO^+(1,3) matrix
% first O(1,3):
% A^{-1}^T V^T eta V A^{-1} = eta
% V^T eta V eta = A^T eta A eta
temp = V.'*eta*V*eta;

% degenerate Sigma^2 => block-diagonal A
% (n-block, entirely in spatial or temporal part:) polar decomposition => A = O(n) * sqrtm(temp)
% (n-block spanning both temporal and spatial parts:) fuck
d = diag(Sigma2);
d = d -[d(2:4); -1];
d = find(abs(d) > tol); % block boundaries
d = d - [0; d(1:end-1)]; % block structure
if (d(1) ~= 1)
    disp('unhandled case, mixed spatial and temporal parts')
    d
    xxx = 1:d(1);
    temp(xxx,xxx)
end

A = sqrtm(temp); % O(n) gauge freedom left
% signs
if (V(1,1) < 0)
  A = A*T; % make it orthochronous
end
if (det(V/A) < 0) % fixme xor
  A = A*eta; % make it proper
end
V = V/A;

% A has still SO^+(1,3)-allowed signature freedom, i.e. two -1s can be multiplied into the space part
Sigma = sqrt(Sigma2); % the square ignores the signature, we still have to restore it.
end


function [v,d] = eigsort(a)
% sorted eig (descending real part)
  [v,d] = eig(a);
  %[v,d] = jordan(a);
  dd = diag(d);
  [~, I] = sort(real(dd), 'descend');
  %d = diag(dd(I));
  d = d(I,I);
  v = v(:,I);
end
