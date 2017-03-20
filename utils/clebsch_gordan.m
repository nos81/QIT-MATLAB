function [ret, JM] = clebsch_gordan(j1, j2, m1, m2, J, M)
% CLEBSCH_GORDAN Clebsch-Gordan coefficients.
%  c = clebsch_gordan(j1, j2, m1, m2, J, M)
%
%  Returns the Clebsch-Gordan coefficient <j1, m1; j2, m2|J, M> used in angular momentum coupling.
%  If only j1 and j2 are given, returns the corresponding unitary basis transform matrix. 

% Ville Bergholm 2017


if nargin == 6
    % return one coefficient
    ret = cg(j1, j2, m1, m2, J, M);
elseif nargin == 2
    % return the entire basis transform matrix
    ret = [];
    JM = [];  % basis ordering
    for J = abs(j1-j2):j1+j2
        temp = [];
        JM = [JM, [J*ones(1,2*J+1); J:-1:-J]];
        for M = -J:J
            for m1 = -j1:j1,
                for m2 = -j2:j2
                    temp((j1-m1)*(2*j2+1) +j2-m2 +1, J-M+1) = cg(j1,j2,m1,m2,J,M);
                end
            end
        end
        ret = [ret, temp];
    end
else
    error('Wrong input syntax.')
end
end


function ret = cg(j1, j2, m1, m2, J, M)
  ret = (-1)^(j1-j2+M) * sqrt(2*J+1) * wigner3j(j1, j2, J, m1, m2, -M);
end
