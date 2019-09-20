function [mBestImage, mBestMembership, bestObjVal, totalIterNum] =...
    binaryMembershipBMAlgor(mAdj, k, convEpsilon, runNum, sUpdateFunc, sDistanceFunc,...
        sUpdateMembershipFunc, sImageInit, sMemInit, varargin)
%
% Matrix factorisation based algorithms for blockmodelling.
%
% Input:
% mAdj          - Adjacency matrix (n x n).
% k             - Number of positions to infer.
% convEpsilon   - Epsilon value to denote convergence.
% runNum        - Number of runs to perform.  
% sUpdateFunc           - Name of the image matrix update function to use.
% sDistanceFunc         - Name of the distance function between A and blockmodel approximation.
% sUpdateMembershipFunc - Name of the membership update function to use.
% sUpdateMembershipFunc - Name of the image intialisation method to use.
% sUpdateMembershipFunc - Name of the membership intialisation method to use.
%
% Output:
% mBestImage        - Image matrix (k x k).
% mBestMembership   - Membership matrix (n x k).
% bestObjVal        - Objective value obtained.
% totalIterNum      - Total number of iterations, over all runs.
%
%
% If using this code, please kindly cite the following paper:
% J. Chan, W. Liu, C. Leckie, J. Bailey and K. Ramamohanarao. 
% "Discovering Latent Blockmodels in Sparse and Noisy Graphs using Non-Negative Matrix Factorisation."
% In Proceedings of 22nd ACM International Conference on Information and Knowledge Management, October 2013.
%
% @author: Jeffrey Chan, 2013
%

    % parse arguments
    inParser = inputParser;
    inParser.KeepUnmatched = true;
    
    defaultSigmoidWeight = size(mAdj,1)^2 / k^2;
    defaultSigmoidMean = 0.5;

    % add parameters and set default valuse
    addParameter(inParser, 'runIterationNum', 100);
    addParameter(inParser, 'sigmoidWeight', defaultSigmoidWeight);
    addParameter(inParser, 'sigmoidMean', defaultSigmoidMean);    

    parse(inParser, varargin{:});
    
    runIterNum = inParser.Results.runIterationNum;
    sigmoidRegWeight = inParser.Results.sigmoidWeight;
    sigmoidMean = inParser.Results.sigmoidMean;

    
    % determine update function for image matrix
    switch sUpdateFunc
        case 'projGradDescM'
            fUpdateFunc = ImageProjGradDesc;
        case 'coordDescM'
            fUpdateFunc = ImageCoordDesc;
        % coordindate descent with 0-1 regularisation on M
        case 'coordDescBMM'
            fUpdateFunc = ImageCoordDescBM;    
        otherwise
            error('binaryMembershipBMAlgor:sUpdateFunc', '%s Invalid update function option', sUpdateFunc);
    end


    % determine distance function
    switch sDistanceFunc
        case 'euclidean'
            fDistanceFunc = EucTriFact;
        case 'euclideanAdj'
            expectedDensity = (sum(sum(mAdj)) / size(mAdj,1)^2);
            mAdjWeight = mAdj - expectedDensity;
            fDistanceFunc = EucTriFactAdj(mAdjWeight);
        case 'bmEuclidean'
            fImageReg = SigmoidImageReg(sigmoidMean);
            fDistanceFunc = EucTriFactBM(fImageReg, sigmoidRegWeight);  
        case 'bmEuclideanAdj'
            expectedDensity = (sum(sum(mAdj)) / size(mAdj,1)^2);
            mAdjWeight = mAdj - expectedDensity;

            fImageReg = SigmoidImageReg(sigmoidMean);
            fImageReg.adjustIdealShift(expectedDensity * ones(k, k));

            fDistanceFunc = EucTriFactBMAdj(mAdjWeight, fImageReg, sigmoidRegWeight);
        otherwise
            error('binaryMembershipBMAlgor:sDistanceFunc', 'Invalid distance function option');
    end


    % determine membership update rules
    switch sUpdateMembershipFunc
        case 'hardCIncr'
            fUpdateMembershipFunc = HardMembershipInc;
        case 'hardCIncrAdjEqual'
            mAdjWeight = mAdj - (sum(sum(mAdj)) / size(mAdj,1)^2);
            fUpdateMembershipFunc = HardMembershipIncAdj(mAdjWeight);    
        otherwise
            error('binaryMembershipBMAlgor:sUpdateMembershipFunc', 'Invalid membership update function option');
    end

    
    % image initialisation
    switch sImageInit
        case 'randomInit'
            fImageInitFunc = InitRandomImage(false);
        case 'randomHardInit'
            fImageInitFunc = InitRandomImage(true);
        otherwise 
            error('binaryMembershipMatFactMult:sImageInit', 'Invalid image initialisation option');
    end % end of switch    
    
    
    % membership initialisation
    switch sMemInit
        case 'randomInit'
            fMemInitFunc = InitRandomMem(false);
        case 'randomHardInit'
            fMemInitFunc = InitRandomMem(true);
        otherwise 
            error('binaryMembershipMatFactMult:sImageInit', 'Invalid image initialisation option');
    end % end of switch  


    totalIterNum = 0;

    % initial run
    if (runNum > 0)
        fprintf('run %d\n', 1);
        
        [mBestImage, mBestMembership, bestObjVal, iterNum] =...
            singleRun(mAdj, k, fUpdateFunc, fImageInitFunc, fMemInitFunc, fDistanceFunc,...
                fUpdateMembershipFunc, convEpsilon, runIterNum, varargin{:});

        totalIterNum = totalIterNum + iterNum;
    else
        error('binaryMembershipBMAlgor:runNum', 'runNum is 0 or less, meaning algorithm was not executed for even one run');
    end


    % loop through remaining number of runs
    for r = 2 : runNum
        fprintf('run %d\n', r);
    
%         [mCurrImage, mCurrMembership, currObjVal, currTotalVal, currMatApproxVal, iterNum] =...        
        [mCurrImage, mCurrMembership, currObjVal, iterNum] =...
            singleRun(mAdj, k, fUpdateFunc, fImageInitFunc, fMemInitFunc, fDistanceFunc,...
                fUpdateMembershipFunc, convEpsilon, runIterNum, varargin{:});
            
        totalIterNum = totalIterNum + iterNum;
        
        % see if this is best solution and update as appropriate
        if (currObjVal < bestObjVal)
            bestObjVal = currObjVal;
            mBestImage = mCurrImage;
            mBestMembership = mCurrMembership;
        end        
    end % end of outer for

end % end of function binarykMeansBMAlgor()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mImage, mMembership, currObjVal, iteration] =...
    singleRun(mAdj, k, fUpdateFunc, fImageInitFunc, fMemInitFunc, fDistanceFunc,...
        fUpdateMembershipFunc, convEpsilon, runIterNum, varargin)    
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
        

    % initial random membership matrix
    if ~bInitMem
        mMembership = fMemInitFunc.initMembership(mAdj, k);
    else
        mMembership = mInitMem;
    end

    % initialise mImage
    if ~bInitImage
        mImage = fImageInitFunc.initImage(mAdj, mMembership);    
    else
        mImage = mInitImage;
    end

        
    mMembership = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
    
    mImage = fUpdateFunc.updateImage(mAdj, mImage, mMembership, fDistanceFunc);    
    
    iteration = 1;
    
    % compute objective
    currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
        
    % make sure initial oldCurrVal is larger than the stopping criteria
    oldObjVal = 3 * currObjVal;
    
    % loop till convergence  
    iteration = iteration + 1;
    if runIterNum <= 0
        while (oldObjVal - currObjVal > convEpsilon * oldObjVal)
            
            mMembership = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
        
            mImage = fUpdateFunc.updateImage(mAdj, mImage, mMembership, fDistanceFunc);
            
            % update objective
            oldObjVal = currObjVal;
            currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
            
            iteration = iteration + 1;
        end % end of while
    else
        while iteration <= runIterNum

            mMembership = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
            
            mImage = fUpdateFunc.updateImage(mAdj, mImage, mMembership, fDistanceFunc);    
    
            % update objective
            currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
   
            iteration = iteration + 1;
        end
        
        
    end
    


end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


