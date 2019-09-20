function [ccmResults] =  testDataWithEdgePrediction(sBaseDir, csFilenames, sOutDir, sOutPrefix,...
    vPosRange, numOfFolds, convEpsilon, algorRunNum, runsPerFoldSet, bDiscretiseMembership)
%
% Loads and performs k-fold cross validation of edge prediction.
%

    bDonotDiscretise = false;
    initPosNum = 1;
    
    
%     ccAlgors = {  
%         {'matApprox', 'bolongM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise}
%         };    
    
    ccAlgors = {  
        {'matApprox', 'hardM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
        {'matApprox', 'bolongM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
        {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
        {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCIncr', bDonotDiscretise},...
        {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCIncr', bDonotDiscretise},...        
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise},...
        {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise},...
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjDeg', 'hardCIncrAdjDeg', bDonotDiscretise},...
        {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjDeg', 'hardCIncrAdjDeg', bDonotDiscretise},...        
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise},...
        {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise},...
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'hardCIncrAdjDeg', bDonotDiscretise},...
        {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'hardCIncrAdjDeg', bDonotDiscretise},...
        {'matApprox', 'bolongM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', bDiscretiseMembership},...
        {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', bDiscretiseMembership},...
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCGradDescAdjEqual', bDiscretiseMembership},...
        {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCCoordDescAdjEqual', bDiscretiseMembership},...
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCGradDescAdjDeg', bDiscretiseMembership}        {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCCoordDescAdjDeg', bDiscretiseMembership},...        
        {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', bDiscretiseMembership},...
        };


%     ccAlgors = {...     
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
%         };



    % load each file
    for i = 1 : length(csFilenames)
        sFilename = fullfile(sBaseDir, csFilenames{i});
        % need to determine if weighted or not
        % no weighted atm
        bWeighted = false;
        bSparse = true;
        bAddOneToIndex = false;
        bHeader = true;
        % pass in some random number, with bHeader = true, the sizes would be
        % read from the header
        vNum = 0;        
        [mGraph, mDim, nDim] = loadMatlabGraph(sFilename, vNum, vNum, bWeighted, bSparse, bAddOneToIndex, bHeader);
        assert(mDim == nDim);
        vNum = mDim;
        
        [vSrc, vTar] = find(mGraph);
        % construct k-folds
        vFoldIndices = crossvalind('Kfold', length(vSrc), numOfFolds);

       
        ccmResults = cell(1, length(ccAlgors));
        resultsNum = 4;
        for c = 1 : length(ccAlgors)
            ccmResults{c} = cell(1, vPosRange(2) - vPosRange(1) + 1);  
            for p = 1 : vPosRange(2) - vPosRange(1) + 1
                ccmResults{c}{p} =  zeros(runsPerFoldSet * numOfFolds, resultsNum);
            end
        end        
        
    
        
        assert(length(vPosRange) == 2);
        assert(vPosRange(1) <= vPosRange(2));
        
        % perform the validation
        for f = 1 : runsPerFoldSet
            for k = 1 : numOfFolds        
        
                vTest = (vFoldIndices == k);
%             vTrain = ~vTest;
                % construct blockmodel on vTrain set of edges
                mTestGraph = sparse(vSrc(vTest), vTar(vTest), ones(length(vSrc(vTest)),1), vNum, vNum);
                mTrainGraph = mGraph - mTestGraph;            
            
                for p = vPosRange(1) : vPosRange(2)
        
                    % Test each algorithm
                    for c = 1 : length(ccAlgors)
                        fprintf('Doing %s-%s-%s.\n', ccAlgors{c}{2}, ccAlgors{c}{6}, ccAlgors{c}{7});
                        [mTrainImage, mTrainMembership, ~, ~, ~, ~, runningTime, totalIterNum] = runBMAlgor(mTrainGraph,...
                            ccAlgors{c}{1}, ccAlgors{c}{2}, p, ccAlgors{c}{4}, ccAlgors{c}{5}, ccAlgors{c}{6}, ccAlgors{c}{7}, ccAlgors{c}{8});
                        % reconstruct the density matrix
                        mTrainImageDensity = getImageDensity(mTrainMembership, mTrainGraph);
%                         plotSingleBlockmodel(mTrainGraph, mTrainMembership, 'whatever');
                        % make predictions about the edges
                        mPredGraph = edgePrediction(mTrainMembership, mTrainImage, ccAlgors{c}{6});
                        % evaluate predictive performance
                        predSse = evalPrediction(mTestGraph, mTrainMembership, mTrainImageDensity);
                
                        % compute ROC
                        [vX,vY, ~, auc] = perfcurve(mTestGraph(:), mPredGraph(:), 1);
                        
%                         plot(vX, vY);

                        % store results
                        ccmResults{c}{p-vPosRange(1)+1}((f-1)*runsPerFoldSet + k,:) = [auc, predSse, runningTime, totalIterNum];
                    end
                end % end of loop through number of positions
    
            end % end of loop through folds
        end % end of loop through fold runs
        
        % construct out file
        sOutFile = fullfile(sOutDir, strcat(sOutPrefix, '_', csFilenames{i}, '.results.csv'));
        fOut = fopen(sOutFile, 'a');
        
        % compute the mean and std. dev.
        for c = 1 : length(ccAlgors)
            for p = 1 : vPosRange(2) - vPosRange(1) + 1
                vMean = mean(ccmResults{c}{p}, 1);
                vStd = std(ccmResults{c}{p});
                % write results to file
                fprintf(fOut, '%s, %s, %s, %d, %.2f, %.2f, %.2f, %.2f, %2.f, %.2f, %.2f, %.2f\n',...
                    ccAlgors{c}{2}, ccAlgors{c}{6}, ccAlgors{c}{7}, vPosRange(1) + p - 1, vMean, vStd);
            end
        end        
        
        fclose(fOut);
    end % end of loop through datasets
    

end % end of function


function mPredGraph =  edgePrediction(mTrainMembership, mTrainImage, sObjective)
%
% Makes a prediction for each edge in mTestGraph.
%
    mPredGraph = mTrainMembership * mTrainImage * mTrainMembership';

end % end of function



function [mDensity] = getImageDensity(mTrainMembership, mTrainGraph)
%
% Obtain the density of the blocks in the train graph.
%
    posNum = size(mTrainMembership,2);
    
    mDensity = zeros(posNum, posNum);
    cVertices = cell(1, posNum);
    
    for p = 1 : posNum
        [vR] = find(mTrainMembership(:,p));
        cVertices{p} = vR;
    end
    
    for r = 1 : posNum
        for c = 1 : posNum
            totalNum = length(cVertices{r}) * length(cVertices{c});
            if totalNum > 0
                mDensity(r,c) = length(find(mTrainGraph(cVertices{r}, cVertices{c}))) / totalNum;
            end
        end
    end

end



function sse = evalPrediction(mTestGraph, mTrainPosMembership, mTrainImage)
%
% Evaluate the prediction accuracy of the blockmodel obtained from training
% folds to the test graph built from the test fold.
%
    % sum of errors
    mApprox = mTestGraph - mTrainPosMembership * mTrainImage * mTrainPosMembership';
    sse = trace(mApprox * mApprox');
end % end of function