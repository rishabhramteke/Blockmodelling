function [mBestImage, mBestMembership, bestObjVal, bestMatApproxVal, totalIterNum] = binaryMembershipNewman(mAdj, k, runNum)
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
%
% @author: Jeffrey Chan, 2013
%


fDistanceFunc = EucTriFactBM;
fMatApproxDisFunc = EucTriFact;

totalIterNum = 0;

% initial run
if (runNum > 0)
    disp(sprintf('run %d', 1));
    [mBestImage, mBestMembership, bestObjVal, bestMatApproxVal, iterNum] = singleRun(mAdj, k, fDistanceFunc, fMatApproxDisFunc);
    totalIterNum = totalIterNum + iterNum;
else
    error('binaryMembershipBMAlgor:runNum', 'runNum is 0 or less, meaning algorithm was not executed for even one run');
end


% loop through remaining number of runs
for r = 2 : runNum
    disp(sprintf('run %d', r));
    
    [mCurrImage, mCurrMembership, currObjVal, bestMatApproxVal, iterNum] = singleRun(mAdj, k, fDistanceFunc, fMatApproxDisFunc);
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


function [mImage, mMembership, currObjVal, currMatApproxVal, iterNum] = singleRun(mAdj, k, fDistanceFunc, fMatApproxDisFunc)
%
% A single run of the algorithm.
%

% MODIFICATION: We should update mImage first, given random initial
% membership, as we might get convergence problems where the mImage forces
% one or more positions to lose all its members


    % initial random solution
    initialise();
    

    % run KL

    
    

    % transpose the solution
    mMembership = mMembership';
    currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
    currMatApproxVal = fMatApproxDisFunc.distance(mAdj, mImage, mMembership);
    
end % end of function



function [] = initialise()
    %
    % Generate an initial solution and compute the degrees of each generated
    % position.
    %


end


