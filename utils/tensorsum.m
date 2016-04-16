function res = tensorsum(varargin)
% TENSORSUM  Like kron but adding instead of multiplying.
    
% Ville Bergholm 2008-2016


%c = log(kron(exp(a), exp(b))); % a perverted way of doing it, the exp overflows...

res = 0;
for k = 1:nargin
    temp = [];
    for a=1:size(res, 1)
        row = [];
        for b=1:size(res, 2)
            foo = res(a,b) +varargin{k};
            row = [row, foo];
        end
        temp = [temp; row];
    end
    res = temp;
end
