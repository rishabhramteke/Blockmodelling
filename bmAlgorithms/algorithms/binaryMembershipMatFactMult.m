function [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, totalIterNum] =...
    binaryMembershipMatFactMult(mAdj, k, runNum, sUpdateFunc, sDistanceFunc,...
        sUpdateMembershipFunc, sImageInit, sMemInit, cfValidMeasFunc, bDiscretiseMembership, bColNormalise, varargin)
%
% Input:
% mAdj      - Adjacency matrix (n x n).
% k         - Number of positions to infer.
% convEpsilon   - Epsilon value to denote convergence.
% runNum    - Number of runs to perform.  Each run we start with a random
%             initial mImage.
% sUpdateFunc - The name of the image matrix update function to use.
% sDistanceFunc - Name of the distance function between A and blockmodel approximation.
% sUpdateMembershipFunc - Name of the membership update function to use.
% cfValidMeasFunc - A cell of function objects that measures the obtained
%                   positioned closeness to the ideal solution over the iterations.
% bDiscretiseMembership - If soft membership update method, whether to
% discretise final memberships or not.
% bColNormalsie - For ding and bolong membership update, whether to normalise
% columns to unit length after each iteration.
% bInitMembership - Whether to randomly initialise the membership as well.

% varargin - Extra parameters that are specific to each algorithm.
%
% Output:
% mImage    - Image matrix (k x k).
% mMembership   - Membership matrix (n x k).
% totalIterNum   - Total number of iterations, over all runs.
%
% Computes the hard approximation of A of blockmodel approximate mMem *
% mImage * mMem'.
% Implements the multiplicative approaches of Bo Long and Chris Ding.
%
%
% @author: Jeffrey Chan, 2014
%

    % parse arguments
    inParser = inputParser;

    % add parameters and set default valuse
    addParameter(inParser, 'runIterationNum', 400);
    % runIterNum - Number of iterations to run.  If ==0, then use convEpsilon for
    %               stopping criteria, otherwise run algorithm for X number of iterations.
    addParameter(inParser, 'rowNormaliseAdhoc', false);
    addParameter(inParser, 'l1MemReg', 0);
    addParameter(inParser, 'l2MemReg', 0);
       
    
    
    parse(inParser, varargin{:});
    
    runIterNum = inParser.Results.runIterationNum;
    bRowNormAdHoc = inParser.Results.rowNormaliseAdhoc;
    l1MemReg = inParser.Results.l1MemReg;
    l2MemReg = inParser.Results.l2MemReg;
    



    % determine update function for image matrix
    switch sUpdateFunc
        case 'bolongHardMemM'
            fUpdateFunc = @updateMultBoLongImageHardMem;
%             fInitFunc = InitUpdateImage(@updateMultBoLongImageImageHardMem);
        case 'bolongSoftMemM'
            fUpdateFunc = @updateMultBoLongImageSoftMem;
%             fInitFunc = InitRandomImage(false);
        case 'dingM'
            fUpdateFunc = @updateMultDingImage;
%             fInitFunc = InitRandomImage(false);
        otherwise
            error('binaryMembershipMatFactMult:sUpdateFunc', 'Invalid update function option');
    end % end of switch

   
    % determine distance function
    switch sDistanceFunc
        case 'euclidean'
            fDistanceFunc = EucTriFact;
            fKKTResidualFunc = EucKKTResidual(EucTriFact);
        otherwise
            error('binaryMembershipMatFactMult:sDistanceFunc', 'Invalid distance function option');
    end


    % determine membership update rules
    switch sUpdateMembershipFunc
        case 'softCMultDing'
            fUpdateMembershipFunc = SoftMembershipMultDing(bColNormalise, bRowNormAdHoc, l1MemReg, l2MemReg);
        case 'softCMultLong'
            fUpdateMembershipFunc = SoftMembershipMultLong(bColNormalise, bRowNormAdHoc, l1MemReg, l2MemReg);
        otherwise
            error('binaryMembershipMatFactMults:sUpdateMembershipFunc', 'Invalid membership update function option');
    end
    
    
    switch sImageInit
        case 'randomInit'
            fImageInitFunc = InitRandomImage(false);
        case 'randomHardInit'
            fImageInitFunc = InitRandomImage(true);
        case 'kmeans'
            
        case 'svd'
            fImageInitFunc = InitSvdImage();
        
        otherwise 
            error('binaryMembershipMatFactMult:sImageInit', 'Invalid image initialisation option');
    end % end of switch    
    
    
    
    switch sMemInit
        case 'randomInit'
            fMemInitFunc = InitRandomMem(false);
        case 'randomHardInit'
            fMemInitFunc = InitRandomMem(true);
        case 'kmeans'
            
        case 'svd'
            fMemInitFunc = InitSvdMem(true);
        case 'randomAcol'
            fMemInitFunc = InitRandomAcolMem();
        otherwise 
            error('binaryMembershipMatFactMult:sImageInit', 'Invalid image initialisation option');
    end % end of switch      

    totalIterNum = 0;
    
    % iteration tracking stats
    cvObjVal = cell(1, runNum);
    cvKKTResidual = cell(1, runNum);
    ccvGroundComparison = cell(1, runNum);
    ccmImage = cell(1, runNum);
    ccmMembership = cell(1, runNum);


    % initial run
    if (runNum > 0)
        fprintf('run %d', 1);
        [mBestImage, mBestMembership, bestObjVal, vObjVal, vKKTResidual, cvComparison, cmImage, cmMembership, iterNum] =...
            singleRun(mAdj, k, fUpdateFunc, fImageInitFunc, fMemInitFunc, fDistanceFunc, fKKTResidualFunc, fUpdateMembershipFunc,...
                bDiscretiseMembership, bColNormalise, runIterNum, cfValidMeasFunc);
            
        totalIterNum = totalIterNum + iterNum;            
            
        cvObjVal{1} = vObjVal;
        cvKKTResidual{1} = vKKTResidual;
        ccvGroundComparison{1} = cvComparison;
        ccmImage{1} = cmImage;
        ccmMembership{1} = cmMembership;
    else
        error('binaryMembershipMatFactMults:runNum', 'runNum is 0 or less, meaning algorithm was not executed for even one run');
    end


    % loop through remaining number of runs
    for r = 2 : runNum
        fprintf('run %d', r);
    
        [mCurrImage, mCurrMembership, currObjVal, vObjVal, vKKTResidual, cvComparison, cmImage, cmMembership, iterNum] =...
            singleRun(mAdj, k, fUpdateFunc, fImageInitFunc, fMemInitFunc, fDistanceFunc, fKKTResidualFunc, fUpdateMembershipFunc,...
                bDiscretiseMembership, bColNormalise, runIterNum, cfValidMeasFunc);
            
        totalIterNum = totalIterNum + iterNum;            
            
        cvObjVal{r} = vObjVal;
        cvKKTResidual{r} = vKKTResidual;
        ccvGroundComparison{r} = cvComparison;   
        ccmImage{r} = cmImage;
        ccmMembership{r} = cmMembership;        
    
        % see if this is best solution and update as appropriate
        if (currObjVal < bestObjVal)
            bestObjVal = currObjVal;
            mBestImage = mCurrImage;
            mBestMembership = mCurrMembership;
        end
    end % end of outer for


end % end of function binaryMembershipMatFactMult()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mImage, mMembership, currObjVal, vObjVal, vKKTResidual, cvComparison, cmImage, cmMembership, iteration] =...
    singleRun(mAdj, k, fUpdateFunc, fImageInitFunc, fMemInitFunc, fDistanceFunc, fKKTResidualFunc,...
        fUpdateMembershipFunc, bDiscretiseMembership, bColNormalise, runIterNum, cfValidMeasFunc)
    %
    % A single run of the algorithm.
    %
    % bDiscretiseMembership - Whether to discretise the membership matrix after
    % the optimal is found.  This is relevant if the position membership update
    % approach is a soft/probablistic one.
    %
    
    
    % store the objective values and kkt residual values after each iteration
    vObjVal = -1 * ones(1, runIterNum);
    vKKTResidual = -1 * ones(1, runIterNum);
    cvComparison = cell(1, length(cfValidMeasFunc));
    for c = 1 : length(cfValidMeasFunc)
        cvComparison{c} = -1 * ones(1, runIterNum);
    end
    cmImage = cell(1, runIterNum);
    cmMembership = cell(1, runIterNum);
    

    % initial random membership matrix
    mMembership = fMemInitFunc.initMembership(mAdj, k);

    % initialise mImage
    mImage = fImageInitFunc.initImage(mAdj, mMembership);
    
    if bColNormalise
        % the mImage is updated from transferring the column normalisation
        % to mMembership.
        [mMembership, mImage] = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
    else
        % we do not update mImage
        mMembership = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
    end
        
    
    iteration = 1;
    
    % compute objective
    currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
    vObjVal(iteration) = currObjVal;
    vKKTResidual(iteration) = fKKTResidualFunc.residual(mAdj, mImage, mMembership);
    for c = 1 : length(cfValidMeasFunc)
        cvComparison{c}(iteration) = cfValidMeasFunc{c}.compare(mMembership);
    end
    cmImage{iteration} = mImage;
    cmMembership{iteration} = mMembership;
    
    % loop for a set number of iterations
    iteration = iteration + 1;
    while iteration <= runIterNum
        if mod(iteration,10) == 0
            fprintf('iteration %d\n', iteration);
        end
        % update mMembership
        mMembership = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
            
        % update mImage
        mImage = fUpdateFunc(mAdj, mImage, mMembership, fDistanceFunc);
    
        % update objective
        currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
        vObjVal(iteration) = currObjVal;
        vKKTResidual(iteration) = fKKTResidualFunc.residual(mAdj, mImage, mMembership);
        for c = 1 : length(cfValidMeasFunc)
            cvComparison{c}(iteration) = cfValidMeasFunc{c}.compare(mMembership);
        end        
        cmImage{iteration} = mImage;
        cmMembership{iteration} = mMembership;        
            
        iteration = iteration + 1;
   end
    
    % see if we need to discretise the results
    if bDiscretiseMembership
       mMembership = discretise(mMembership); 
       % compute the new distances
       currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
       vObjVal(end) = currObjVal;
       vKKTResidual(end) = fKKTResidualFunc.residual(mAdj, mImage, mMembership);
       cmMembership{end} = mMembership;
    end

end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


