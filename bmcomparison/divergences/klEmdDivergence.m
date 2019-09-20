function [dist] = klEmdDivergence(vPmf1, vPmf2, vWeightDistValues)
%
% Computes the EMD-KL divergence between the two weight distributions.
%
% @author: Jeffrey Chan
%

    % ensure vWeightDistValues is a column vector
    if size(vWeightDistValues,2) > size(vWeightDistValues,1)
        vWeightDistValues = vWeightDistValues';
    end
    assert(size(vPmf1,1) >= size(vPmf1,2));
    assert(size(vPmf2,1) >= size(vPmf2,2));
    
    
    % compute cost matrix
    mQ = squareform(pdist(vWeightDistValues, 'cityblock'));
    vF = mQ(:)';
        
    % equality constraints: row (col) sums must equal sum of probabilites of
    % first (second) distribution/frequency vectors
    mAeq1 = repmat(diag(ones(1, length(vWeightDistValues))), 1, length(vWeightDistValues));
    mAeq2 = zeros(length(vWeightDistValues), length(vWeightDistValues)*length(vWeightDistValues));
    for u = 1 : length(vWeightDistValues)
        mAeq2(u, ((u-1) * length(vWeightDistValues) + 1) : (u * length(vWeightDistValues))) = ones(1, length(vWeightDistValues));
    end
    % combine the equality constraints
    mAeq = cat(1, mAeq1, mAeq2);
    vBeq = cat(1, vPmf1, vPmf2);
    
    % lower and upper bound of matching weights [0,1]
    vLb = zeros(length(vWeightDistValues) * length(vWeightDistValues), 1);
    vUb = ones(length(vWeightDistValues) * length(vWeightDistValues), 1);

    [vMatching, dist] = linprog(vF, [], [], mAeq, vBeq, vLb, vUb);
end