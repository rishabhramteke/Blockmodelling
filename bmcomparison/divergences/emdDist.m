function [dist] = emdDist(vWeightDist1, vWeightDist2, vWeightDistValues)
%
% Computes the EMD distance between the two weight distributions.
%
% If using this code, please kindly cite the following paper:
% J. Chan, X. Nguyen, W. Liu, C. Leckie, J. Bailey, K. Ramamohanarao and J. Pei. 
% "Structure-aware Distance Measures for Comparing Clusterings in Graphs". 
% In Proceedings of the 18th Pacific-Asia Conference on Knowledge Discovery and Data Mining, May 2014. 
%
% @author: Jeffrey Chan
%

    % ensure vWeightDistValues is a column vector
    if size(vWeightDistValues,2) > size(vWeightDistValues,1)
        vWeightDistValues = vWeightDistValues';
    end
    assert(size(vWeightDist1,1) >= size(vWeightDist1,2));
    assert(size(vWeightDist2,1) >= size(vWeightDist2,2));
    
    
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
    vBeq = cat(1, vWeightDist1, vWeightDist2);
    
    % lower and upper bound of matching weights [0,1]
    vLb = zeros(length(vWeightDistValues) * length(vWeightDistValues), 1);
    vUb = ones(length(vWeightDistValues) * length(vWeightDistValues), 1);

    options = optimoptions('linprog', 'Display', 'off');
    [vMatching, dist] = linprog(vF, [], [], mAeq, vBeq, vLb, vUb, [], options);
end