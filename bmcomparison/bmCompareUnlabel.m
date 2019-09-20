function [dist, mW, vPosWeights1, vPosWeights2] = bmCompareUnlabel(mMembership1, mMembership2, mImageMat1, mImageMat2, bPosDef, sMatchType)
%
% Align and compute the inter image distance based on the image densities.
%
% bPosDef - boolean, whether we wish to use SVD approximation to force the
% hessian to be positive definite.
% sMatchType - string, specifying which type of matching/comparison to be
% performed.  One of {relative, absolute, nosize}
%

% compute the relative weights of each position in each of the
% blockmodels/co-clusterings
% TODO: assuming hard partitionoing
vPosWeights1 = zeros(size(mImageMat1,1), 1);
vPosWeights2 = zeros(size(mImageMat2,1), 1);
for v = 1 : size(mMembership1,1)
   vPos = find(mMembership1(v,:) > 0);
   vPosWeights1(vPos) = vPosWeights1(vPos) + 1;
end
for v = 1 : size(mMembership2,1)
   vPos = find(mMembership2(v,:) > 0);
   vPosWeights2(vPos) = vPosWeights2(vPos) + 1;
end

switch sMatchType
    case 'relative'
        % normalise by the number of elements in each blockmodel/co-clustering
        vPosWeights1 = vPosWeights1 / size(mMembership1,1);
        vPosWeights2 = vPosWeights2 / size(mMembership2,1);
    case 'absolute'
        % do not need to normalise
    case 'nosize'
        % we assign weights that are equal and equal to the larger sized
        % of the two sets of positions
        % This will force best matchings and any unmatched we can easily
        % detect with total weightings of 0.
        newSize = max(size(mImageMat1,1), size(mImageMat2,1));
        vPosWeight1(:,:) = 1 / newSize;
        vPosWeight2(:,:) = 1 / newSize;
    otherwise
        warning('Unknown sMatchType = %s specified', sMatchType);
        return;
end


% compute the inter position density distances
mInterPosDist = computeInterPosDen(mImageMat1, mImageMat2);

% compute the distance
[dist, mW] = matchingDist(mInterPosDist, vPosWeights1, vPosWeights2, bPosDef, sMatchType);


end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mInterPosDist] = computeInterPosDen(mImageMat1, mImageMat2)
%
% Compute the inter position density distances
%

mInterPosDist = zeros(size(mImageMat1, 1), size(mImageMat1, 2), size(mImageMat2, 1), size(mImageMat2, 2));
for r1 = 1 : size(mImageMat1, 1)
    for c1 = 1 : size(mImageMat1, 2)
        for r2 = 1 :size(mImageMat2, 1)
            for c2 = 1 :size(mImageMat2, 2)
                mInterPosDist(r1,c1,r2,c2) = absDenDist(mImageMat1(r1,c1), mImageMat2(r2,c2));
            end
        end
    end
end

end % end of function


function [dist] = absDenDist(density1, density2)

dist = abs(density1 - density2);

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [dist, vW] = matchingDist(mInterPosDist, vPosWeights1, vPosWeights2, bPosDef, sMatchType)
%
% mInterPosDist - the distances between the positions
% vPosWeight1 - the proportional weights for bm1 (must sum to 1)
% vPosWeight2 - the proportional weights for bm2 (must sum to 1)

assert(size(mInterPosDist,1) == size(vPosWeights1,1));
assert(size(mInterPosDist,2) == size(vPosWeights2,1));
assert(sum(vPosWeights1) == 1);
assert(sum(vPosWeights2) == 1);


k1 = size(vPosWeights1,1);
k2 = size(vPosWeights2,1);
k1k2 = k1 * k2;

%
% write the objective
%

% quadratic terms (Hessian)
mH = zeros(k1k2, k1k2);
for r1 = 1 : k1
    for c1 = 1 : k1
        for r2 = 1 : k2
            for c2 = 1: k2
                mH((r1 -1) * k2 + r2, (c1-1) * k2 + c2) = mInterPosDist(r1,c1,r2,c2);
            end
        end
    end
end
% make mH symmetric first
mH = (mH + mH') / 2;

if bPosDef
    % need to compute the eigenvalues of mH and ensure they are all positive
    % to ensure mH is positive (semi)definite
    [mEigenVecs, mEigenVals] = eig(mH);
    % ensure all mEigenVals are (semi)positive
    mEigenVals(mEigenVals < 0) = eps;
    % reconstruct matrix
    mH = mEigenVecs * mEigenVals * mEigenVecs';
end




% (linear terms (we have none, so set to zero vector)
%mF = zeros(k1k2,1);


%
% write the linear constraints
%

mAeq = [];
vBeq = [];
switch sMatchType
    case 'absolute'
        [mAExtra, vBExtra] = conEquConstraintsRel(vPosWeights1, vPosWeights2, k1, k2, k1k2);
        mA = [mA ; mAExtra];
        vB = [vB ; vBExtra];
    otherwise
        % equality constraints
        % \sum W = alpha or beta's
        [mAeq, vBeq] = conEquConstraintsRel(vPosWeights1, vPosWeights2, k1, k2, k1k2);
end


% intial values (k1 x k2) (solve simultaneous equations AX = B)
mXInt = mAeq \ vBeq;

% set the options
options = optimset('Algorithm','active-set');
fObjective = @(x)(alignmentObj(x, mH));
[vW, dist, exitFlag] = fmincon(fObjective, mXInt, mA, vB, mAeq, vBeq, [], [], [], options);

% vW = quadprog(mH, mF, mA, vB, mAeq, vBeq);
% dist = vW' * mH * vW;

end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mAeq, vBeq] = conEquConstraintsRel(vPosWeights1, vPosWeights2)
%
% constructs the equility constraints (to do with position sizes) based on
% the relative size idea.
%

% k1 = size(vPosWeights1,1);
% k2 = size(vPosWeights2,1);
% k1k2 = k1 * k2;


mAeq = zeros(k1+k2, k1k2);
vBeq = zeros(k1+k2,1);
% do alpha's first
for i = 1 : k1
    for j = 1 : k2
       mAeq(i,(i-1)*k2 + j) = 1; 
    end
    vBeq(i) = vPosWeights1(i);
end

% now do beta's
for j = 1 : k2
    for i = 1 : k1
       mAeq(k1+j, (i-1)*k2 + j) = 1; 
    end
    vBeq(k1+j) = vPosWeights2(j);
end
% project to lower space
AeqRank = rank(mAeq);
if (AeqRank < size(mAeq,1))
    [U,S,V] = svd(mAeq);
    mAeq = S([1:AeqRank],[1:AeqRank]) * V(:,[1:AeqRank])';
    vBeq = U' * vBeq;
    vBeq = vBeq([1:AeqRank]);
end

end % end of function



function [mAeq, vBeq] = conEquConstraintsAbs(vPosWeights1, vPosWeights2)
%
% constructs the equility constraints (to do with position sizes) based on
% the absolute size idea.
%

% k1 = size(vPosWeights1,1);
% k2 = size(vPosWeights2,1);
% k1k2 = k1 * k2;


mAeq = zeros(k1+k2, k1k2);
vBeq = zeros(k1+k2,1);
% do alpha's first
for i = 1 : k1
    for j = 1 : k2
       mAeq(i,(i-1)*k2 + j) = 1; 
    end
    vBeq(i) = vPosWeights1(i);
end

% now do beta's
for j = 1 : k2
    for i = 1 : k1
       mAeq(k1+j, (i-1)*k2 + j) = 1; 
    end
    vBeq(k1+j) = vPosWeights2(j);
end
% project to lower space
AeqRank = rank(mAeq);
if (AeqRank < size(mAeq,1))
    [U,S,V] = svd(mAeq);
    mAeq = S([1:AeqRank],[1:AeqRank]) * V(:,[1:AeqRank])';
    vBeq = U' * vBeq;
    vBeq = vBeq([1:AeqRank]);
end

end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f = alignmentObj(vW, mH)
f = vW' * mH * vW;
end


% function f = objective(X, mInterPosDist)
% 
% f = 0.0;
% 
% for r1 = size(mInterPosDist, 1)
%     for c1 = size(mInterPosDist, 2)
%         for r2 = size(mInterPosDist, 1)
%             for c2 = size(mInterPosDist, 2)
%                 f = f + X(r1,r2) * X(c1,c2) * mInterPosDist(r1,c1,r2,c2);
%             end
%         end
%     end
% end


% end % end of function
