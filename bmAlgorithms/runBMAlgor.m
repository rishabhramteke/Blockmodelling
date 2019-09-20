% function [mBestImage, mBestMembership, bestObjVal, bestMatApproxVal, idealEucObjVal, idealManhObjVal, runningTime, totalIterNum] = runBMAlgor(mAdj,...
%     sObj, sAlgor, posNum, convEpsilon, runNum, sDist, sPosAlgor, bDiscretiseMembership, bColNormalise, varargin)
function [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, runningTime, totalIterNum] =...
    runBMAlgor(mAdj, sObj, sAlgor, posNum, convEpsilon, runNum, sDist, sPosAlgor, sImageInit, sMemInit, bDiscretiseMembership, bColNormalise, fImageDistanceFunc, cfValidMeasFunc, varargin)

%
% Runs the specified blockmodelling algorithm.
%
% sObj -        The objective to use.
% sAlgor -      The optimisation (image) algorithm.
% posNum -      Number of positions.
% convEpsilon - Stopping threshold (epsilon).
% runNum -      Number of runs to execute.
% sDist -       The distance function.
% sPosAlgor -   Membership optimisation algorithm.
% bDiscretiseMembership - Whether to discretise the membership matrix after
% the optimal is found.  This is relevant if the position membership update
% approach is a soft/probablistic one.
% varargin -    Additional parameters that are passed in
%               binaryMembershipBMAlgor.
%
%
% @author: Jeffrey Chan, 2013
%

    cvObjVal = {}; cvKKTResidual = {}; ccvGroundComparison = {}; ccmImage = {}; ccmMembership = {};

    t = tic();
    
    switch sObj
        case 'matApprox'
            % two arguments in varargin
            %         assert(length(varargin) == 2);
            [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, totalIterNum] =...
                binaryMembershipBMAlgor(mAdj, posNum, convEpsilon, runNum, sAlgor, sDist, fImageDistanceFunc, sPosAlgor, sImageInit, sMemInit, cfValidMeasFunc, bDiscretiseMembership, bColNormalise, varargin{:});
%         case 'matFact'
%             [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, totalIterNum] =...
%                 binaryMembershipMatFactMult(mAdj, posNum, runNum, sAlgor, sDist, sPosAlgor, sImageInit, sMemInit, cfValidMeasFunc, bDiscretiseMembership, bColNormalise, varargin{:});         
        case 'hardSa'
            [mBestImage, mBestMembership, bestObjVal, bestMatApproxVal, totalIterNum] = binaryMembershipSA(mAdj, posNum, runNum, sDist);
        case 'reichardt'
            [mBestImage, mBestMembership, bestObjVal, totalIterNum] = binaryMembershipReichardtSA(mAdj, posNum, runNum, sDist, varargin{:});
        case 'irm'
            [mBestImage, mBestMembership, bestObjVal, totalIterNum] = binaryMembershipIRM(mAdj, posNum, runNum);
        case 'bnmtf'
            [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, totalIterNum] =...
                binaryMembershipBNMTF(mAdj, posNum, runNum, fImageDistanceFunc, cfValidMeasFunc, bDiscretiseMembership, varargin{:});
        case 'bkn'
            % poisson generative model of Ball, Karrer and Newman
            [mBestImage, mBestMembership, bestObjVal, totalIterNum] =...
                binaryMembershipBBK(mAdj, posNum, runNum, convEpsilon, varargin{:});      
        case 'mmsbm'
            % mixed membership stochastic blockmodel
            [mBestImage, mBestMembership, bestObjVal, totalIterNum] =...
                binaryMembershipMMSB(mAdj, posNum, runNum, convEpsilon);                     
        otherwise
            error('runBMAlgor:sObj', 'Invalide objective function name');     
    end % end of switch

    runningTime = toc(t);


    % compare mBestImage with ideal
%     [mIdeal, idealEucObjVal] = computeEucIdealDist(mBestImage);
%     idealManhObjVal = computeManhIdealDist(mBestImage);


end % end of function


function [mIdeal, idealObjVal] = computeEucIdealDist(mImage)
%
% Computes the eucldiean distance between M_ideal and mImage.
%

mIdeal = zeros(size(mImage,1), size(mImage,2));
for r = 1 : size(mImage,1)
    for c = 1 : size(mImage,2)
       if (1-mImage(r,c)) <= mImage(r,c)
           mIdeal(r,c) = 1;
       end
    end
end

mDiff = mIdeal - mImage;
idealObjVal = trace(mDiff * mDiff');

end % end of function


function [idealObjVal] = computeManhIdealDist(mImage)
%
% Computes the manhanttan distance between M_ideal and mImage.
%

mIdeal = zeros(size(mImage,1), size(mImage,2));
for r = 1 : size(mImage,1)
    for c = 1 : size(mImage,2)
       if (1-mImage(r,c)) <= mImage(r,c)
           mIdeal(r,c) = 1;
       end
    end
end

idealObjVal = sum(sum(abs(mIdeal - mImage)));

end % end of function


