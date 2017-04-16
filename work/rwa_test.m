function L = rwa_test(H, D, baths)
% RWA_TEST  Test the validity of the rotating wave approximation in
% the derivation of the weak-coupling Born-Markov Lindblad equation.
%
%  In the weak-coupling case the Lindblad equation is obtained
%  after invoking a rotating wave approximation which discards all rotating
%  terms, all assumed to rotate so rapidly that this is justified.
%
%  The dH differences must not be too low for the RWA to hold.
%  The smallest dH difference determines \tau_S, the intrinsic time scale of
%  the system, as defined in Breuer&Petruccione, chapter 3.3.1.
%  For the RWA to work, \tau_S must be much shorter than \tau_R, the
%  relaxation time of the system.
%
%  Probably more accurately, the rotation speed of each term must be
%  much more than its norm, or maybe the norm of its restriction to
%  the space of hermitian matrices?
%
%  B can also be a cell vector of baths, in which case D has to be
%  a cell vector of the corresponding interaction operators.
%  The baths are assumed independent and uncorrelated,
%  i.e. \Gamma_{ab} = \delta_{ab} * \Gamma_{a}.

% Ville Bergholm 2017


[dH, X, G] = markov.components(H, D, baths);
w = dH/2/pi
n_baths = size(X,2);

for n=1:n_baths
    for j=1:length(dH)
        normA(j,n) = norm(X{j,n});
    end
end
normA

% norms of the jump ops
tol_norm_A = 1e-10;

% Swapping w and w' gives the (counter-rotating) h.c. term, therefore
% it's enough to sum over w < w'. The diagonal w = w' contains the
% nonrotating terms.

% TODO some w'-w rates may be accidentally identical, in this case
% the operators should be summed before taking the norm.

diag = {};
antidiag = {};
antidiag_dw = [];
nondiag = {};
dw = []; % rotation rate: w'-w

% NOTE we assume uncorrelated baths, hence just a single sum!
for n=1:n_baths
    A = X(:,n);  % ops for the n'th bath
    % j: w
    for j=1:length(dH)
        j
        w = dH(j);
        if normA(j,n) < tol_norm_A
            continue
        end
        % nonrotating diagonal term, (w,w)
        diag{end+1} = op(G(n,j,1), G(n,j,1), A{j}, A{j});

        if w ~= 0
            % nonrotating diagonal term, (-w,-w)
            diag{end+1} = op(G(n,j,2), G(n,j,2), A{j}', A{j}');
            
            % antidiagonal (-w, w) term
            antidiag_dw(end+1) = 2*w;
            antidiag{end+1} = op(G(n,j,2), G(n,j,1), A{j}', A{j});
        end
        % k: w'
        for k=j+1:length(dH)
            z = dH(k);
            if normA(k,n) < tol_norm_A
                continue
            end
            % negative dH vals: vary both separately
            % we always get a pair of nondiagonal ops rotating at the same rate: (w,z) ^= (-z,-w)
            dw(end+1) = z-w;
            temp = op(G(n,j,1), G(n,k,1), A{j}, A{k}) +op(G(n,k,2), G(n,j,2), A{k}', A{j}');
            nondiag{end+1} = temp;
            if w ~= 0
                % another nondiag pair: (-w,z) ^= (-z,w)
                dw(end+1) = z+w;
                temp = op(G(n,j,2), G(n,k,1), A{j}', A{k}) +op(G(n,k,2), G(n,j,1), A{k}', A{j});
                nondiag{end+1} = temp;
            end
        end
    end
end

% sum up the diagonal, nonrotating component
L = 0;
for k=1:length(diag)
    L = L + diag{k};
end

% TODO instead of norm(), some special norm on the hermitian restriction?
figure;
hold on
xlabel('|w''-w|')
ylabel('|dL|')

dw = dw/2/pi
antidiag_dw = antidiag_dw/2/pi

plot(0, norm(L), 'ko')
plot(antidiag_dw, map_cellvec(antidiag, @norm), 'bo')
plot(dw, map_cellvec(nondiag, @norm), 'ro')
grid on
legend('diag', 'antidiag', 'nondiag')
end


function ret = op(Ga, Gb, A, B)
  % build the rotating lindblad term
  T = B'*A;
  % TEST discard extra levels, AFTER computing T to get Lamb shift correct?
  %keep = [1:3, 5:7, 9:11];
  %T = T(keep, keep);
  %A = A(keep, keep);
  %B = B(keep, keep);
  
  ret = (Ga +conj(Gb)) * lrmul(A,B') -Ga * lmul(T) -conj(Gb) * rmul(T);
  % TODO project op into hermitian subspace (it should be HP anyway)
end


function ret = map_cellvec(v, f)
% Apply a function to every element of a cell vector.
    ret = [];
    for k=1:length(v)
        ret(k) = f(v{k});
    end
end
