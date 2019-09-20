function [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, runningTime, totalIterNum] =...
    binaryMembershipHyrbid(mAdj, sObj1, sObj2, sAlgor1, sAlgor2, posNum, convEpsilon1, convEpsilon2, runNum, sDist1, sDist2, sPosAlgor1, sPosAlgor2,...
        sImageInit1, sMemInit1, bDiscretiseMembership1, bDiscretiseMembership2, bColNormalise1, bColNormalise2, fImageDistanceFunc, cfValidMeasFunc, cVarargin1, cVarargin2)
    %
    % Hybrid algorithm.
    %
    % If set number of iterations, should set up separately in cVarargin's.
    %
    % @author: Jeffrey Chan, 2014
    %
    
    
    
    % algorithm run number is only one
    algorRunNum = 1;
    
    % first algorithm
    [mBestImage1, mBestMembership1, bestObjVal1, cvObjVal1, cvImageDis1, cvKKTResidual1,...
        ccvGroundComparison1, ccmImage1, ccmMembership1, runningTime1, totalIterNum1] =...
            runBMAlgor(mAdj, sObj1, sAlgor1, posNum, convEpsilon1, algorRunNum, sDist1,...
                sPosAlgor1, sImageInit1, sMemInit1, bDiscretiseMembership1, bColNormalise1, fImageDistanceFunc, cfValidMeasFunc, cVarargin1{:});                    

    imageSum = sum(sum(mBestImage1));
    mDiag = diag(repmat(imageSum, 1, posNum));
    mInvDiag = diag(repmat(1/imageSum,1, posNum));
    mBestImage1 = mInvDiag * mBestImage1;
    mBestMembership1 = mBestMembership1 * mDiag;
            
    % second algorithm
    % use the sImageInit1 and sMemInit1, as place holders as the membership will
    % be set to the best one from previous run
    [mBestImage, mBestMembership, bestObjVal2, cvObjVal2, cvImageDis2, cvKKTResidual2,...
        ccvGroundComparison2, ccmImage2, ccmMembership2, runningTime2, totalIterNum2] =...
            runBMAlgor(mAdj, sObj2, sAlgor2, posNum, convEpsilon2, algorRunNum, sDist2,...
                sPosAlgor2, sImageInit1, sMemInit1, bDiscretiseMembership2, bColNormalise2, fImageDistanceFunc, cfValidMeasFunc,...
                    'initialMem', mBestMembership1, 'initialImage', mBestImage1, cVarargin2{:});       
                

    % TODO: check these objective values are the same or at least same scaling
%     if bestObjVal1 < bestObjVal2
%         mBestImage = mBestImage1;
%         mBestMembership = mBestMembership1;
%         bestObjVal = bestObjVal1;
%     else
%         bestObjVal = bestObjVal2;
%     end
    bestObjVal = bestObjVal2;
                
    % concatenate the collected statistics
    cvObjVal = cell(1, algorRunNum);
    for i = 1 : algorRunNum
        cvObjVal{i} = cat(2, cvObjVal1{i}, cvObjVal2{i});
    end
    

    cvImageDis = cell(1, algorRunNum);
    for i = 1 : algorRunNum
        cvImageDis{i} = cat(2, cvImageDis1{i}, cvImageDis2{i});
    end
    

    cvKKTResidual = cell(1, algorRunNum);
    for i = 1 : algorRunNum
        cvKKTResidual{i} = cat(2, cvKKTResidual1{i}, cvKKTResidual2{i});
    end
    

    ccvGroundComparison = cell(1, algorRunNum);
    for i = 1 : algorRunNum
        ccvGroundComparison{i} = cell(1, size(ccvGroundComparison1{i},2));
        for j = 1 : size(ccvGroundComparison1{i}, 2)
            ccvGroundComparison{i}{j} = cat(2, ccvGroundComparison1{i}{j}, ccvGroundComparison2{i}{j});
        end
    end
    

    ccmImage = cell(1, algorRunNum);
    for i = 1 : algorRunNum
        ccmImage{i} = cat(2, ccmImage1{i}, ccmImage2{i});
    end    

    ccmMembership = cell(1, algorRunNum);
    for i = 1 : algorRunNum
        ccmMembership{i} = cat(2, ccmMembership1{i}, ccmMembership2{i});
    end      
    
    runningTime = runningTime1 + runningTime2;
    totalIterNum = totalIterNum1 + totalIterNum2;

end % end of function

