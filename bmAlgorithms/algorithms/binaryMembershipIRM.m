function [mBestImage, mBestMembership, bestObjVal, totalIterNum] = binaryMembershipIRM(mAdj, k, runNum)
%
% Input:
% mAdj      - Adjacency matrix (n x n).
% k         - Number of positions to infer.
% convEpsilon   - Epsilon value to denote convergence.
% runNum    - Number of runs to perform.  Each run we start with a random
%             initial mImage.
% sUpdateFunc - The name of the update function to use.
%
% Output:
% mImage    - Image matrix (k x k).
% mMembership   - Membership matrix (n x k).
%
% Computes the hard approximation of A of blockmodel approximate mMem *
% mImage * mMem'.
% Uses Joerg Reichardt's null model formulation.
%


    fDistanceFunc = EucTriFactBM;

    totalIterNum = 0;

    % initial run
    if (runNum > 0)
        fprintf('run %d', 1);
        [mBestImage, mBestMembership, bestObjVal, iterNum] = singleRun(mAdj, k, fDistanceFunc);
        totalIterNum = totalIterNum + iterNum;
    else
        error('binaryMembershipBMAlgor:runNum', 'runNum is 0 or less, meaning algorithm was not executed for even one run');
    end


    % loop through remaining number of runs
    for r = 2 : runNum
        fprintf('run %d', r);

        [mCurrImage, mCurrMembership, currObjVal, iterNum] = singleRun(mAdj, k, fDistanceFunc);
        totalIterNum = totalIterNum + iterNum;

        % see if this is best solution and update as appropriate
        if (currObjVal < bestObjVal)
            bestObjVal = currObjVal;
            mBestImage = mCurrImage;
            mBestMembership = mCurrMembership;
        end
    end % end of outer for




end % end of function binarykMeansBMAlgor()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555


function [mImage, mMembership, currObjVal, iterNum] = singleRun(mAdj, k, fDistanceFunc)
%
% A single run of the algorithm.
%


    % parameters
    opts = struct();
    % only one relation, so a one by one cell
    A = cell(1, 1);
    W = cell(1,1);
    A{1} = mAdj;
    W{1} = zeros(size(mAdj,1), size(mAdj,2));
    
    % run IRM inference algorithm
    [L,cpu_time,mMembership,mImage,sample,West] = IRMUnipartite(A, W, k, opts);
    
    iterNum = size(L,1);
    

    % transpose the solution
    mMembership = mMembership';
    currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
    
end % end of function




