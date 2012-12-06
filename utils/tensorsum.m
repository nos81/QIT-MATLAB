function c = tensorsum(a, b)
% TENSORSUM  Like kron but adding instead of multiplying.
    
% Ville Bergholm 2008-2010

    
%c = log(kron(exp(a), exp(b))); % a perverted way of doing it, the exp overflows...
c = [];
for k=1:length(a)
  c = [c, a(k)+b];
end
