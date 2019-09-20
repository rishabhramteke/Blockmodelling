function [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, totalIterNum] =...
    binaryMembershipBNMTF(mAdj, k, runNum, fImageDistanceFunc, cfValidMeasFunc, bDiscretiseMembership, varargin)
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
% @author: Jeffrey Chan, 2013
%


    % parse arguments
    inParser = inputParser;
    inParser.KeepUnmatched = true;

    % add parameters and set default valuse
    addParameter(inParser, 'lossOption', 0);
    addParameter(inParser, 'sparseOption', 0);
    addParameter(inParser, 'regWeight', 0);    
    addParameter(inParser, 'maxRank', k);
    addParameter(inParser, 'runIterationNum', 100);
    % default is to collect per iteration stats
    addParameter(inParser, 'collectPerIterStats', true);    


    parse(inParser, varargin{:});
    
    bLossOption = inParser.Results.lossOption;
    bSparseOption = inParser.Results.sparseOption;
    regWeight = inParser.Results.regWeight;
    maxRank = inParser.Results.maxRank;    
    maxIter = inParser.Results.runIterationNum;   
    bCollectPerIterStats = inParser.Results.collectPerIterStats;  

    


    % iteration tracking stats
    cvObjVal = cell(1, runNum);
    cvImageDis = cell(1, runNum);
    cvKKTResidual = cell(1, runNum);
    ccvGroundComparison = cell(1, runNum);
    ccmImage = cell(1, runNum);
    ccmMembership = cell(1, runNum);


    totalIterNum = 0;
    
    fprintf('%s\n', 'bnmtfM');

    % initial run
    if (runNum > 0)
        fprintf('run %d\n', 1);
        [mBestImage, mBestMembership, bestObjVal, vObjVal, vImageDis, vKKTResidual, cvComparison, cmImage, cmMembership, iterNum] =...
            singleRun(mAdj, fImageDistanceFunc, cfValidMeasFunc, bDiscretiseMembership, bLossOption, bSparseOption, regWeight, maxRank, maxIter, varargin{:});
    
        totalIterNum = totalIterNum + iterNum;
        
        if bCollectPerIterStats
            cvObjVal{1} = vObjVal;
            cvImageDis{1} = vImageDis;
            cvKKTResidual{1} = vKKTResidual;
            ccvGroundComparison{1} = cvComparison;
            ccmImage{1} = cmImage;
            ccmMembership{1} = cmMembership;        
        end
    else
        error('binaryMembershipBNMTF:runNum', 'runNum is 0 or less, meaning algorithm was not executed for even one run');
    end


    % loop through remaining number of runs
    for r = 2 : runNum
        fprintf('run %d\n', r);
    
        [mCurrImage, mCurrMembership, currObjVal, vObjVal, vImageDis, vKKTResidual, cvComparison, cmImage, cmMembership, iterNum] =...
            singleRun(mAdj, fImageDistanceFunc, cfValidMeasFunc, bDiscretiseMembership, bLossOption, bSparseOption, regWeight, maxRank, maxIter, varargin{:});
    
        totalIterNum = totalIterNum + iterNum;
        
        if bCollectPerIterStats
            cvObjVal{r} = vObjVal;
            cvImageDis{r} = vImageDis;
            cvKKTResidual{r} = vKKTResidual;
            ccvGroundComparison{r} = cvComparison;   
            ccmImage{r} = cmImage;
            ccmMembership{r} = cmMembership;            
        end
    
        % see if this is best solution and update as appropriate
        if (currObjVal < bestObjVal)
            bestObjVal = currObjVal;
            mBestImage = mCurrImage;
            mBestMembership = mCurrMembership;
        end
    end % end of outer for

end % end of function binarykMeansBMAlgor()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555


function [mImage, mMembership, currObjVal, vObjVal, vImageDis, vKKTResidual, cvComparison, cmImage, cmMembership, iterNum] =...
    singleRun(mAdj, fImageDistanceFunc, cfValidMeasFunc, bDiscretiseMembership, bLossOption, bSparseOption, regWeight, maxRank, maxIter, varargin)
    %
    % A single run of the algorithm.
    %
    
    % parse arguments
%     inParser = inputParser;   
%     inParser.KeepUnmatched = true;
% 
% 
%     parse(inParser, varargin{:});
     

    
    
    % run the BNMTF algorithm
    % we do not specify the initial matrices
    % default will run for a max iteration of 100
    [structModel, currObjVal, vObjVal, vImageDis, vKKTResidual, cvComparison, cmImage, cmMembership] =...
        BNMTF(mAdj, bLossOption, bSparseOption, regWeight, fImageDistanceFunc, cfValidMeasFunc, bDiscretiseMembership, maxRank, maxIter, varargin{:});          

    
    mImage = structModel.B;
    mMembership = structModel.U;

    iterNum = maxIter;
    
end % end of function



