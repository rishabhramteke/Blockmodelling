function [mBestImage, mBestMembership, bestObjVal, totalIterNum] =...
    binaryMembershipBBK(mAdj, k, runNum, convEpsilon, varargin)
%
% 
% Implement the EM algorithm of Ball, Karrer and Newman in the paper "An
% efficient and principled method for detecting communities in networks".
%
% Input:
% mAdj      - Adjacency matrix (n x n).
% k         - Number of positions to infer.
% runNum    - Number of runs to perform.  Each run we start with a random
%             initial mImage.
% convEpsilon
%
% Output:
% mBestImage        - Image matrix (k x k).
% mBestMembership   - Membership matrix (n x k).
% bestObjVal        - 
% bestMatApproxVal  - 
% totalIterNum      - 
%


    fDistanceFunc = EucTriFactBM;
    
    % initialisation for membership
    fMemInitFunc = InitRandomMem(false);
    fImageInitFunc = InitRandomImage(false);
    

    totalIterNum = 0;

    % initial run
    if (runNum > 0)
        fprintf('run %d\n', 1);
        [mBestImage, mBestMembership, bestObjVal, iterNum] = singleRun(mAdj, k, convEpsilon, fDistanceFunc, fMemInitFunc, fImageInitFunc, varargin{:});
        totalIterNum = totalIterNum + iterNum;
    else
        error('binaryMembershipBBK:runNum', 'runNum is 0 or less, meaning algorithm was not executed for even one run');
    end


    % loop through remaining number of runs
    for r = 2 : runNum
        fprintf('run %d\n', r);
    
        [mCurrImage, mCurrMembership, currObjVal, iterNum] = singleRun(mAdj, k, convEpsilon, fDistanceFunc, fMemInitFunc, fImageInitFunc, varargin{:});
        totalIterNum = totalIterNum + iterNum;
    
        % see if this is best solution and update as appropriate
        if (currObjVal < bestObjVal)
            bestObjVal = currObjVal;
            mBestImage = mCurrImage;
            mBestMembership = mCurrMembership;
        end
    end % end of outer for

end % end of function binaryMembershipBBK()


function [mImage, mMembership, currObjVal, iterNum] =...
    singleRun(mAdj, k, convEpsilon, fDistanceFunc, fMemInitFunc, fImageInitFunc, varargin)
%
% A single run of the algorithm.
%

    % parse arguments
    inParser = inputParser;
    inParser.KeepUnmatched = true;

    % add parameters and set default valuse
    % initial membership and image matrices
    addParameter(inParser, 'initialMem', NaN);
    addParameter(inParser, 'initialImage', NaN);
    
    parse(inParser, varargin{:});
        
    mInitMem = inParser.Results.initialMem;
    mInitImage = inParser.Results.initialImage;

    bInitMem = false;
    bInitImage = false;
    if ~isnan(mInitMem)
        bInitMem = true;
    end
    if ~isnan(mInitImage)
        bInitImage = true;
    end
    
    
    % run the EM algorithm
    iterNum = 1;
    
    % initialise membership and image
    if ~bInitMem
        mMembership = fMemInitFunc.initMembership(mAdj, k);
    else
        mMembership = mInitMem;
    end
    if ~bInitImage
        mImage = fImageInitFunc.initImage(mAdj, mMembership);    
    else
        mImage = mInitImage;
    end    
 
    
    prevLikelihood = logLikelihood(mAdj, mMembership, mImage);
    
    % compute initial first loop
    mQ = zeros(size(mAdj,1), size(mAdj,1), k, k);
    [mMembership, mImage, mQ] = oneLoop(mAdj, mMembership, mImage, mQ);
    
    % new likelihood
    currLikelihood = logLikelihood(mAdj, mMembership, mImage);
    
    while currLikelihood - prevLikelihood > convEpsilon
        iterNum = iterNum + 1;
        prevLikelihood = currLikelihood;
        [mMembership, mImage, mQ] = oneLoop(mAdj, mMembership, mImage, mQ);
        currLikelihood = logLikelihood(mAdj, mMembership, mImage);
    end    
    
    
    % calculate how far
    currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
    
end % end of function



function [likelihood] = logLikelihood(mAdj, mMembership, mImage)
%
% Compute the log likelihood.
%

    mPosSum = mMembership * mImage * mMembership';
    likelihood = sum(sum(mAdj .* log(mPosSum))) - sum(sum(mPosSum));

end % end of function


function [mMembership, mImage, mQ] = oneLoop(mAdj, mMembership, mImage, mQ)
    
    k = size(mMembership,2);
    % update mQ
    mPosSum = mMembership * mImage * mMembership';
    for i = 1 : size(mAdj,1)
        for j = 1 : size(mAdj,1)
            for r = 1 : k
                for s = 1 : k
                    mQ(i,j,r,s) = mMembership(i,r) * mImage(r,s) * mMembership(j,s) / mPosSum(i,j);
                end
            end
        end
    end
    
    % update membership
    mNomin = zeros(size(mAdj,1),k);
    for r = 1 : k
        for i = 1 : size(mAdj,1)
            for j = 1 : size(mAdj,1)
                for s = 1 : k
                    mNomin(i,r) = mNomin(i,r) + mAdj(i,j) * mQ(i,j,r,s);
                end
            end       
        end
    end
    
    for r = 1 : k
        mMembership(:,r) = mNomin(:,r) / sum(mNomin(:,r));
    end
    
    % update mImage
    for r = 1 : k
        for s = 1 : k
            mImage(r,s) = sum(sum(mAdj .* mQ(:,:,r,s)));
        end
    end


end