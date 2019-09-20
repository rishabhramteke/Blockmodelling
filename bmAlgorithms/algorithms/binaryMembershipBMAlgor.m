function [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, totalIterNum] =...
    binaryMembershipBMAlgor(mAdj, k, convEpsilon, runNum, sUpdateFunc, sDistanceFunc, fImageDistanceFunc,...
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
% Matrix factorisation based algorithms for blockmodelling.
%
%
% @author: Jeffrey Chan, 2013
%

    % parse arguments
    inParser = inputParser;
    inParser.KeepUnmatched = true;
    
    defaultSigmoidWeight = 0.1 * size(mAdj,1)^2 / k^2;
    defaultSigmoidMean = 0.5;

    % add parameters and set default valuse
    addParameter(inParser, 'runIterationNum', 100);
    addParameter(inParser, 'sigmoidWeight', defaultSigmoidWeight);
    addParameter(inParser, 'sigmoidMean', defaultSigmoidMean);    
    addParameter(inParser, 'rowNormaliseAdhoc', false);
    % regularisation parameters for multipllicative algorithms
    addParameter(inParser, 'l1MemReg', 0);
    addParameter(inParser, 'l2MemReg', 0);    
    addParameter(inParser, 'l1ImageReg', 0);
    addParameter(inParser, 'l2ImageReg', 0);    
    % default is to collect per iteration stats
    addParameter(inParser, 'collectPerIterStats', true);
    addParameter(inParser, 'penaltyWeight', 0.01);
            
            

   
    parse(inParser, varargin{:});
    
    runIterNum = inParser.Results.runIterationNum;
    sigmoidRegWeight = inParser.Results.sigmoidWeight;
    sigmoidMean = inParser.Results.sigmoidMean;
    bRowNormAdHoc = inParser.Results.rowNormaliseAdhoc;
    l1MemReg = inParser.Results.l1MemReg;
    l2MemReg = inParser.Results.l2MemReg;    
    l1ImageReg = inParser.Results.l1ImageReg;
    l2ImageReg = inParser.Results.l2ImageReg;   
    bCollectPerIterStats = inParser.Results.collectPerIterStats;   
    penaltyWeight = inParser.Results.penaltyWeight;        
    
   
    fprintf('%s + %s + %s\n', sUpdateFunc, sUpdateMembershipFunc, sDistanceFunc);

    % determine update function for image matrix
    switch sUpdateFunc
        case 'hardM'
            fUpdateFunc = @updateHardImage;
%             fInitFunc = InitRandomImage;
        case 'kmeansM'
            fUpdateFunc = @updateAvgImage;
%             fInitFunc = @initBootstrapImage;
        case 'bolongHardMemM'
%             fUpdateFunc = @updateMultBoLongImageHardMem;
            fUpdateFunc = ImageMultLongHardMem;
%             fInitFunc = InitUpdateImage(@updateMultBoLongImageImageHardMem);
        case 'bolongSoftMemM'
%             fUpdateFunc = @updateMultBoLongImageSoftMem;
            fUpdateFunc = ImageMultLongSoftMem(l1ImageReg, l2ImageReg);
%             fInitFunc = InitRandomImage(false);
        case 'dingM'
%             fUpdateFunc = @updateMultDingImage;
            fUpdateFunc = ImageMultDing(l1ImageReg, l2ImageReg);
%             fInitFunc = InitRandomImage(false);
        case 'projGradDescM'
%             fUpdateFunc = @updateProjGradDescImage;
            fUpdateFunc = ImageProjGradDesc;
%             fInitFunc = InitRandomImage(false);
        case 'coordDescM'
            fUpdateFunc = ImageCoordDesc;
%             fInitFunc = InitRandomImage(false);
        % coordindate descent with 0-1 regularisation on M
        case 'coordDescBMM'
            fUpdateFunc = ImageCoordDescBM;
%             fInitFunc = InitRandomImage(false);        
        % exact step size (for euclidean problem only)
        case 'coordDescExactM'
%             fUpdateFunc = @updateCoordDescImageExact;
            fUpdateFunc = ImageCoordDescExact(true, l1ImageReg);
%             fInitFunc = InitRandomImage(false);                
        case 'coordDescExactBMM'
%             fUpdateFunc = @updateCoordDescImageExact;
            fUpdateFunc = ImageCoordDescBMExact(l1ImageReg);
%             fInitFunc = InitRandomImage(false);        
        case 'multBMM'
            fUpdateFunc = ImageMultBM(l1ImageReg, l2ImageReg);
        % multiplicative update (with potential normalisation and 0-1
        % regularisation)
        case 'multReparaM'
            fUpdateFunc = ImageMultHardSumConstraintRepara(l1ImageReg, l2ImageReg);
        % langragian constraint of M sum
        case 'multLangM'            
            fUpdateFunc = ImageMultHardSumConstraint(l1ImageReg, l2ImageReg);
        otherwise
            error('binaryMembershipBMAlgor:sUpdateFunc', '%s Invalid update function option', sUpdateFunc);
    end


% determine distance function
switch sDistanceFunc
    case 'euclidean'
        fDistanceFunc = EucTriFact;
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);
%         fMatApproxDisFunc = EucTriFact;
%         fTotalDisFunc = EucTriFact;
    case 'euclideanBinPen'
        fPenalty = BinaryMembershipPenaltyTerm(penaltyWeight);
        fDistanceFunc = EucTriFact('computeExtraTerms', fPenalty);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);
    case 'euclideanAdjEqual'
        expectedDensity = (sum(sum(mAdj)) / size(mAdj,1)^2);
        mAdjWeight = mAdj - expectedDensity;
        fDistanceFunc = EucTriFactAdj(mAdjWeight);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);
%         fDistanceFunc.adjustIdealShift(expectedDensity * ones(k, k));
%         fMatApproxDisFunc = EucTriFact;
%         fTotalDisFunc = fDistanceFunc;
    case 'euclideanAdjEquaBinPen'
        expectedDensity = (sum(sum(mAdj)) / size(mAdj,1)^2);
        mAdjWeight = mAdj - expectedDensity;
        fPenalty = BinaryMembershipPenaltyTerm(penaltyWeight);
        fDistanceFunc = EucTriFactAdj(mAdjWeight, 'computeExtraTerms', fPenalty);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);
    case 'euclideanAdjDeg'
        vInDeg = sum(mAdj,1);
        vOutDeg = sum(mAdj,2);
        % outer product
        edgeNum = sum(sum(mAdj));
        mExpectedDensity = vInDeg * vOutDeg / (edgeNum^2);
        mAdjWeight = mAdj - mExpectedDensity;
        fDistanceFunc = EucTriFactAdjDeg(mAdjWeight);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);
%         fDistanceFunc.adjustIdealShift(mExpectedDensity);
%         fMatApproxDisFunc = EucTriFact;
%         fTotalDisFunc = fDistanceFunc;           
    case 'bmEuclidean'
        fImageReg = SigmoidImageReg(sigmoidMean);
        fDistanceFunc = EucTriFactBM('imageRegFunc', fImageReg, 'imageRegWeight', sigmoidRegWeight);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);
        
%         fMatApproxDisFunc = EucTriFact;
%         fTotalDisFunc = fDistanceFunc;
    case 'bmEuclideanBinPen'
        fImageReg = SigmoidImageReg(sigmoidMean);
        fPenalty = BinaryMembershipPenaltyTerm(penaltyWeight);
        fDistanceFunc = EucTriFactBM('imageRegFunc', fImageReg, 'imageRegWeight', sigmoidRegWeight, 'computeExtraTerms', fPenalty);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);
    case 'bmEuclideanAdjEqual'
        expectedDensity = (sum(sum(mAdj)) / size(mAdj,1)^2);
        mAdjWeight = mAdj - expectedDensity;
        
        fImageReg = SigmoidImageReg(sigmoidMean);
        fImageReg.adjustIdealShift(expectedDensity * ones(k, k));
        fDistanceFunc = EucTriFactBMAdj(mAdjWeight, 'imageRegFunc', fImageReg, 'imageRegWeight', sigmoidRegWeight);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);
%         fDistanceFunc.adjustIdealShift(expectedDensity * ones(k, k));
%         fMatApproxDisFunc = EucTriFact;
%         fTotalDisFunc = fDistanceFunc;
    case 'bmEuclideanAdjEqualBinPen'
        expectedDensity = (sum(sum(mAdj)) / size(mAdj,1)^2);
        mAdjWeight = mAdj - expectedDensity;
        
        fImageReg = SigmoidImageReg(sigmoidMean);
        fImageReg.adjustIdealShift(expectedDensity * ones(k, k));
        
        fPenalty = BinaryMembershipPenaltyTerm(penaltyWeight);
        
        fDistanceFunc = EucTriFactBMAdj(mAdjWeight, 'imageRegFunc', fImageReg, 'imageRegWeight', sigmoidRegWeight, 'computeExtraTerms', fPenalty);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);
%         fDistanceFunc.adjustIdealShift(expectedDensity * ones(k, k));
%         fMatApproxDisFunc = EucTriFact;
%         fTotalDisFunc = fDistanceFunc;
    case 'bmEuclideanAdjDeg'
        vInDeg = sum(mAdj,1);
        vOutDeg = sum(mAdj,2);
        % outer product
        edgeNum = sum(sum(mAdj));
        mExpectedDensity = vInDeg * vOutDeg / (edgeNum^2);
        mAdjWeight = mAdj - mExpectedDensity;
            
        fImageReg = SigmoidImageReg(sigmoidMean);
        fImageReg.adjustIdealShift(mExpectedDensity);

        fDistanceFunc = EucTriFactBMAdjDeg(mAdjWeight, fImageReg, sigmoidRegWeight);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);
      
%         fMatApproxDisFunc = EucTriFact;
%         fTotalDisFunc = fDistanceFunc;      
    case 'bmCubicEuclidean'
        fImageReg = CubicImageReg(sigmoidMean * ones(k,k));
        fDistanceFunc = EucTriFactBM(fImageReg, sigmoidRegWeight);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);        
    case 'bmCubicEuclideanAdjEqual'
        expectedDensity = (sum(sum(mAdj)) / size(mAdj,1)^2);
        mAdjWeight = mAdj - expectedDensity;
        
        fImageReg = CubicImageReg( expectedDensity * ones(k,k) );
%         fImageReg.adjustIdealShift(expectedDensity * ones(k, k));
        fDistanceFunc = EucTriFactBMAdj(mAdjWeight, fImageReg, sigmoidRegWeight);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);
    case 'bmCubicEuclideanAdjDeg'
        vInDeg = sum(mAdj,1);
        vOutDeg = sum(mAdj,2);
        % outer product
        edgeNum = sum(sum(mAdj));softCMultDing
        mExpectedDensity = vInDeg * vOutDeg / (edgeNum^2);
        mAdjWeight = mAdj - mExpectedDensity;
            
        fImageReg = CubicImageReg(mExpectedDensity);
%         fImageReg.adjustIdealShift(mExpectedDensity);

        fDistanceFunc = EucTriFactBMAdjDeg(mAdjWeight, fImageReg, sigmoidRegWeight);
        fKKTResidualFunc = EucKKTResidual(fDistanceFunc);        
    otherwise
        error('binaryMembershipBMAlgor:sDistanceFunc', 'Invalid distance function option');
end




% determine membership update rules
switch sUpdateMembershipFunc
    case 'hardCBatch'
        fUpdateMembershipFunc = HardMembershipBatch;
    case 'hardCIncr'
        fUpdateMembershipFunc = HardMembershipInc;
    case 'hardCIncrAdjEqual'
        mAdjWeight = mAdj - (sum(sum(mAdj)) / size(mAdj,1)^2);
        fUpdateMembershipFunc = HardMembershipIncAdj(mAdjWeight);
    case 'hardCIncrAdjDeg'
        vInDeg = sum(mAdj,1);
        vOutDeg = sum(mAdj,2);
        % outer product
        mAdjWeight = mAdj - vInDeg * vOutDeg / sum(sum(mAdj));
        fUpdateMembershipFunc = HardMembershipIncAdj(mAdjWeight);      
    case 'hardCPenaltyRowConstraint'
        fUpdateMembershipFunc = HardMembershipPenaltyRowConstraint('initPenaltyWeight', penaltyWeight, 'normaliseCol', bColNormalise, 'l1RegWeight',  l1MemReg, 'l2RegWeight', l2MemReg);                
    case 'hardCPenaltyRowConstraintAdjEqual'
        mAdjWeight = mAdj - (sum(sum(mAdj)) / size(mAdj,1)^2);
        fUpdateMembershipFunc = HardMembershipPenaltyRowConstraint('initPenaltyWeight', penaltyWeight, 'normaliseCol', bColNormalise, 'l1RegWeight',  l1MemReg, 'l2RegWeight', l2MemReg,'adjWeight', mAdjWeight);                        
    case 'softCMultDing'
        fUpdateMembershipFunc = SoftMembershipMultDing(bColNormalise, bRowNormAdHoc, l1MemReg, l2MemReg);
    case 'softCMultLong'
        fUpdateMembershipFunc = SoftMembershipMultLong(bColNormalise, bRowNormAdHoc, l1MemReg, l2MemReg);
    case 'softCMultHardRowConstraint'
        fUpdateMembershipFunc = SoftMembershipMultHardRowConstraint('normaliseCol', bColNormalise, 'l1RegWeight',  l1MemReg, 'l2RegWeight', l2MemReg);        
    case 'softCMultHardRowConstraintAdjEqual'
        mAdjWeight = mAdj - (sum(sum(mAdj)) / size(mAdj,1)^2);
        fUpdateMembershipFunc = SoftMembershipMultHardRowConstraint('normaliseCol', bColNormalise, 'l1RegWeight',  l1MemReg, 'l2RegWeight', l2MemReg, 'adjWeight', mAdjWeight);
    case 'softCMultHardRowConstraintRepara'
        fUpdateMembershipFunc = SoftMembershipMultHardRowConstraintRepara(bColNormalise, l1MemReg, l2MemReg);                
    % soft with orthogonal constraints on mMembership
    case 'softCMultOrthogonal'
        fUpdateMembershipFunc = SoftMembershipMultOrthogonal();                
    case 'softCCoordDesc'
        fUpdateMembershipFunc = SoftMembershipCoord;
    case 'softCGradDesc'
        fUpdateMembershipFunc = SoftMembershipGrad(fDistanceFunc);
    case 'softCCoordDescAdjEqual'
        mAdjWeight = mAdj - (sum(sum(mAdj)) / size(mAdj,1)^2);
        fUpdateMembershipFunc = SoftMembershipCoordAdj(mAdjWeight);
    case 'softCCoordDescAdjDeg'
        vInDeg = sum(mAdj,1);
        vOutDeg = sum(mAdj,2);
        % outer product
        mAdjWeight = mAdj - vInDeg * vOutDeg / sum(sum(mAdj));
        fUpdateMembershipFunc = SoftMembershipCoordAdj(mAdjWeight);  
    case 'softCGradDescAdjEqual'
        % altered by appropriate fDistanceFunc (need to ensure that
        % fDistanceFunc is the correct one)
        fUpdateMembershipFunc = SoftMembershipGrad(fDistanceFunc);        
    case 'softCGradDescAdjDeg'
        % altered by appropriate fDistanceFunc (need to ensure that
        % fDistanceFunc is the correct one)
        fUpdateMembershipFunc = SoftMembershipGrad(fDistanceFunc);    
    otherwise
        error('binaryMembershipBMAlgor:sUpdateMembershipFunc', 'Invalid membership update function option');
end

    % image initialisation
    switch sImageInit
        case 'randomInit'
            fImageInitFunc = InitRandomImage(false);
        case 'randomHardInit'
            fImageInitFunc = InitRandomImage(true);
        case 'kmeans'
            fImageInitFunc = InitKmeansImage();
        case 'svd'
            fImageInitFunc = InitSvdImage();
        otherwise 
            error('binaryMembershipMatFactMult:sImageInit', 'Invalid image initialisation option');
    end % end of switch    
    
    
    % membership initialisation
    switch sMemInit
        case 'randomInit'
            fMemInitFunc = InitRandomMem(false);
        case 'randomHardInit'
            fMemInitFunc = InitRandomMem(true);
        case 'kmeans'
            fMemInitFunc = InitKmeansMem(true);
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
    cvImageDis = cell(1, runNum);
    cvKKTResidual = cell(1, runNum);
    ccvGroundComparison = cell(1, runNum);
    ccmImage = cell(1, runNum);
    ccmMembership = cell(1, runNum);    
    
    
    % initial runcfValidMeasFunc
    if (runNum > 0)
        fprintf('run %d\n', 1);
        
        [mBestImage, mBestMembership, bestObjVal, vObjVal, vImageDis, vKKTResidual, cvComparison, cmImage, cmMembership, iterNum] =...
            singleRun(mAdj, k, fUpdateFunc, fImageInitFunc, fMemInitFunc, fDistanceFunc, fImageDistanceFunc, fKKTResidualFunc,...
                fUpdateMembershipFunc, bDiscretiseMembership, bColNormalise, convEpsilon, runIterNum, cfValidMeasFunc, bCollectPerIterStats, varargin{:});

            
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
        error('binaryMembershipBMAlgor:runNum', 'runNum is 0 or less, meaning algorithm was not executed for even one run');
    end


    % loop through remaining number of runs
    for r = 2 : runNum
        fprintf('run %d\n', r);
    
%         [mCurrImage, mCurrMembership, currObjVal, currTotalVal, currMatApproxVal, iterNum] =...        
        [mCurrImage, mCurrMembership, currObjVal, vObjVal, vImageDis, vKKTResidual, cvComparison, cmImage, cmMembership, iterNum] =...
            singleRun(mAdj, k, fUpdateFunc, fImageInitFunc, fMemInitFunc, fDistanceFunc, fImageDistanceFunc, fKKTResidualFunc,...
                fUpdateMembershipFunc, bDiscretiseMembership, bColNormalise, convEpsilon, runIterNum, cfValidMeasFunc, bCollectPerIterStats, varargin{:});
            
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [mImage, mMembership, currObjVal, currTotalVal, currMatApproxVal, iteration] =...
%     singleRun(mAdj, k, convEpsilon, fUpdateFunc, fInitFunc, fDistanceFunc, fMatApproxDisFunc,...
%         fTotalDisFunc, fUpdateMembershipFunc, bDiscretiseMembership, bColNormalise, runIterNum)
function [mImage, mMembership, currObjVal, vObjVal, vImageDis, vKKTResidual, cvComparison, cmImage, cmMembership, iteration] =...
    singleRun(mAdj, k, fUpdateFunc, fImageInitFunc, fMemInitFunc, fDistanceFunc, fImageDistanceFunc, fKKTResidualFunc,...
        fUpdateMembershipFunc, bDiscretiseMembership, bColNormalise, convEpsilon, runIterNum, cfValidMeasFunc, bCollectPerIterStats, varargin)    
    %
    % A single run of the algorithm.
    %
    % bDiscretiseMembership - Whether to discretise the membership matrix after
    % the optimal is found.  This is relevant if the position membership update
    % approach is a soft/probablistic one.
    % bCollectPerIterStats  - Whether to collect iteration by iteration stats.
    %   If false, then vObjVal, vKKTResidual, cmImage, cvComparison,
    %   cmMembership are empty.
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
    

    % store the objective values and kkt residual values after each iteration
    vObjVal = -1 * ones(1, runIterNum);
    vImageDis = -1 * ones(1, runIterNum);
    vKKTResidual = -1 * ones(1, runIterNum);
    cvComparison = cell(1, length(cfValidMeasFunc));
    for c = 1 : length(cfValidMeasFunc)
        cvComparison{c} = -1 * ones(1, runIterNum);
    end
    cmImage = cell(1, runIterNum);
    cmMembership = cell(1, runIterNum);    
    

    
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

    % set zero entries to epsilon
    mMembership = max(eps, mMembership);
    mImage = max(eps, mImage);
        

    if bColNormalise
        % the mImage is updated from transferring the column normalisation
        % to mMembership.
        [mMembership, mImage] = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
    else
        % we do not update mImage
        mMembership = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
    end
    
    mImage = fUpdateFunc.updateImage(mAdj, mImage, mMembership, fDistanceFunc);    
    
    iteration = 1;
    
    % compute objective
    currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
    
    if bCollectPerIterStats
        vObjVal(iteration) = currObjVal;
         vImageDis(iteration) = fImageDistanceFunc.distance(mImage, mMembership);
        vKKTResidual(iteration) = fKKTResidualFunc.residual(mAdj, mImage, mMembership);
        for c = 1 : length(cfValidMeasFunc)
            cvComparison{c}(iteration) = cfValidMeasFunc{c}.compare(mMembership);
        end
        cmImage{iteration} = mImage;
        cmMembership{iteration} = mMembership;
    end
    
%     currMatApproxVal = fMatApproxDisFunc.distance(mAdj, mImage, mMembership);
%     currTotalVal = fTotalDisFunc.distance(mAdj, mImage, mMembership);
    % make sure initial oldCurrVal is larger than the stopping criteria
    oldObjVal = 3 * currObjVal;
    
    % loop till convergence  
    iteration = iteration + 1;
    if runIterNum <= 0
        while (oldObjVal - currObjVal > convEpsilon * oldObjVal)
            fprintf('iteration %d\n', iteration);
            
            
            
            if bColNormalise
                % the mImage is updated from transferring the column normalisation
                % to mMembership.
                [mMembership, mImage] = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
            else
                % we do not update mImage
                mMembership = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
            end            
            
            mImage = fUpdateFunc.updateImage(mAdj, mImage, mMembership, fDistanceFunc);
            
   
            % update objective
            oldObjVal = currObjVal;
            currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);

            if bCollectPerIterStats
                vObjVal(iteration) = currObjVal;
                vImageDis(iteration) = fImageDistanceFunc.distance(mImage);
                vKKTResidual(iteration) = fKKTResidualFunc.residual(mAdj, mImage, mMembership);
                for c = 1 : length(cfValidMeasFunc)
                    cvComparison{c}(iteration) = cfValidMeasFunc{c}.compare(mMembership);
                end        
                cmImage{iteration} = mImage;
                cmMembership{iteration} = mMembership;   
            end
    
            iteration = iteration + 1;
        end % end of while
    else
        while iteration <= runIterNum
%             if mod(iteration,10) == 0
%                 fprintf('iteration %d\n', iteration);
%             end
            
                    
            
            if bColNormalise
                % the mImage is updated from transferring the column normalisation
                % to mMembership.
                [mMembership, mImage] = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
            else
                % we do not update mImage
                mMembership = fUpdateMembershipFunc.updateMembership(mAdj, mImage, mMembership);
            end      
            
            mImage = fUpdateFunc.updateImage(mAdj, mImage, mMembership, fDistanceFunc);    
    

    
            % update objective
            currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);

            if bCollectPerIterStats
                vObjVal(iteration) = currObjVal;
                vImageDis(iteration) = fImageDistanceFunc.distance(mImage, mMembership);
                vKKTResidual(iteration) = fKKTResidualFunc.residual(mAdj, mImage, mMembership);
                for c = 1 : length(cfValidMeasFunc)
                    cvComparison{c}(iteration) = cfValidMeasFunc{c}.compare(mMembership);
                end        
                cmImage{iteration} = mImage;
                cmMembership{iteration} = mMembership;               
            end
    
            iteration = iteration + 1;
        end
        
        
    end
    
    % see if we need to discretise the results
    if bDiscretiseMembership
       mMembership = discretise(mMembership); 
       % compute the new distances
       currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
%        currMatApproxVal = fMatApproxDisFunc.distance(mAdj, mImage, mMembership);
%        currTotalVal = fTotalDisFunc.distance(mAdj, mImage, mMembership);
       vObjVal(end) = currObjVal;
       if bCollectPerIterStats
           vKKTResidual(end) = fKKTResidualFunc.residual(mAdj, mImage, mMembership);
           cmMembership{end} = mMembership;       
       end
    end

end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% 
% function [mImage] = initEmptyImage(mAdj, mMembership)
% %
% % Initialises the image matrix.
% %
% 
%     mImage = [];
% 
% end 



function [mImage] = initBootstrapImage(mAdj, mMembership)
%
% Initialises the image matrix.
%
% Use hard initialisation.
%

    mImage = updateAvgImage(mAdj, [], mMembership);

end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5








function [mImage] = updateHardImage(mAdj, mImage, mMembership, fDistanceFunc)
%
% Updates mImage using EM/kMeans like local optimisation.
% Each individual optimisation is independent of the what values the other
% densities are set to, hence the optimisation can be order-independent.
%

% calls updateAvgImage, which basically computes the density matrix using
% current mMembership.  We then set M(row,col) to 0 F_norm(row,col) <= 0.5
% or 1 otherwise.

mDensity = updateAvgImage(mAdj, mImage, mMembership);
% intialise mImage to all zeros
mImage = zeros(size(mDensity,1), size(mDensity,2));
[vR, vC, vVal] = find(mDensity > 0.5);
% assign all entries with densities greater than 0.5 to 1
mImage(vR, vC) = 1;


end % end of function



function [mImage] = updateAvgImage(mAdj, mImage, mMembership,fDistanceFunc)
%
% Updates mImage using EM/kMeans like local optimisation.
% Each individual optimisation is independent of the what values the other
% densities are set to, hence the optimisation can be order-independent.
%

% divide each column of mMembership by the scalar of its sum.  This is like
% normalising wrt to the number of elements in each position
% When we compute F = C'*A*C, F is kxk, which is number of ones in each
% block.  When we normalise C as above, this gives us a F(row,col) that is
% normalised by |C_row||C_col|, or the density of the block (or average,
% like kMeans).  We set M(row,col) to F_norm(row,col), and this will
% minimise the euclidean difference between the block and approximated
% block.

% normalise mMembership
mNormMembership = mMembership;
vColNorm = sum(mMembership,1);
% make sure all values of vColNum > 0
assert(size(find(vColNorm), 2) == size(mMembership,2));

for c = 1 : size(mMembership, 2)
   vNonZero = find(mMembership(:,c));
   mNormMembership(vNonZero, c) = mNormMembership(vNonZero, c) / vColNorm(c);
end

mImage = mNormMembership' * mAdj * mNormMembership;


end % end of function








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


