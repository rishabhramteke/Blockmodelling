
function [value] = myEntropy(vDistribution)
%
% compute entropy
%


value = 0;
for p = 1 : size(vDistribution,2)
    if vDistribution(p) > 0
        value = value + vDistribution(p) * log2(vDistribution(p));
    end
end

value = -value;

end