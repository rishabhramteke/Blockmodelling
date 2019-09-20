function [stResults, mBestMembership, mBestImage] = evalSingleDataset(cData, varargin)
%
% Evaluate the scaling issues with factorisation algorithms.
%
% @author: Jeffrey Chan, 2014
%
    
%
% Scans the input directories then runs the testing
%
% cData -
% numRum - the number of runs to experiment over ([0,numRum]).  The appropriate
% input run files must exist.
% sOutBaseDir - Base directory of output files.
% sOutPrefix - prefix of output files.
% 
% Options (varargin):
% 
%
%

    % parse arguments
    inParser = inputParser;
    inParser.KeepUnmatched = true;

    % add parameters and set default valuse
    addParameter(inParser, 'discreteMemb', false);
    addParameter(inParser, 'normaliseCol', false);
    addParameter(inParser, 'convergEsp', 0.01);
%     addParameter(inParser, 'posBinNum', 100);
%     addParameter(inParser, 'imageBinNum', 10);
    

    parse(inParser, varargin{:});
    
    bDiscretiseMembership = inParser.Results.discreteMemb;
    bColNorm = inParser.Results.normaliseCol;
    convEpsilon = inParser.Results.convergEsp;
%     posBinNum = inParser.Results.posBinNum;
%     imageBinNum = inParser.Results.imageBinNum;
    
    
    
    % self specify
    algorRunNum = 1;
    initPosNum = 1; %this will be updated with the actual posNum when the ground truths are loaded


%             {'s', {'matApprox', 'multReparaM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0}},...
%             {'s', {'matApprox', 'multReparaM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0, 'l1MemReg', 100, 'l1ImageReg', 10}},...
%             {'s', {'matApprox', 'multLangM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0, 'l1MemReg', 100}},...            
%             {'s', {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.5}},...                    
%                 {'s', {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 400}},...
%                 {'s', {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon+1, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 400, 'rowNormaliseAdhoc', true}},...
%                cAlgor = {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 400};
%                cAlgor = {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclideanAdjEqual', 'softCMultHardRowConstraintAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.0,  'runIterationNum', 100};             
%                 cAlgor = {'matApprox', 'dingM', initPosNum, convEpsilon+1, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 400, 'rowNormaliseAdhoc', true};
%                 cAlgor = {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.0,  'runIterationNum', 100};
                cAlgor = {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraintRepara', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.0,  'runIterationNum', 100};
%                 {'s', {'matApprox', 'multLangM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0}},...
%                 cAlgor = {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'regWeight', 0, 'runIterationNum', 100, 'innerCollectPerIterStats', false};                        
%                 {'s', {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclideanAdjEqual', 'softCMultHardRowConstraintAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.0,  'runIterationNum', 100}},...             
%                 {'s', {'matApprox', 'multLangM', initPosNum, convEpsilon+1, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0, 'l1MemReg', 100}},...                                            
%             {'s', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 10}},... 
%             {'s', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 10}},...
%             {'s', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultHardRowConstraintRepara', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 100}},...              
%             {'h', {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'regWeight', 0, 'runIterationNum', 10},  {'matApprox', 'multLangM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0, 'runIterationNum', 1}},...                        

    

    % initialisations (same for all algorithms)
    sImageInit = 'randomInit';
    sMemInit = 'randomInit';    

    [stResults, mBestMembership, mBestImage] = performRuns(cAlgor, cData, sImageInit, sMemInit, varargin{:});
    
    

    
    % close the files
%     for c = 1 : size(cData,2)
%         fclose(cfOutResults{c});
%     end
    
        
end % end of function



function [stResults, mBestMembership, mBestImage] = performRuns(cAlgors, cData, sImageInit, sMemInit, varargin)


    % parse arguments
%     inParser = inputParser;
%     inParser.KeepUnmatched = true;    
% 
%     % default is to collect per iteration stats
%     addParameter(inParser, 'collectPerIterStats', true);
%     
%     parse(inParser, varargin{:});
    
%     bCollectPerIterStats = inParser.Results.collectPerIterStats;   
    

       
    

    
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
                    error('runAndCompare:sImageInit', 'Invalid image initialisation option');
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
                    error('runAndCompare:sMemInit', 'Invalid membership initialisation option');
    end % end of switch             
    
    



    % PARAMETERS

    

    mAdj = cData{1};
    numPos = cData{2};
    sDataName = cData{3};
%     numPos = size(mVertPosMapAct,2);
                

    % update number of positions in the cccAlgors
    cAlgors{3} = numPos;
  
        
%     % comparison measures between ground truth and obtained positions
%     if bCollectPerIterStats
%         cfValidMeasFunc = {...
%             ClusterComparisonWrapper(@(a,b)(fuzzyComparison(a,b,'rand')), mVertPosMapAct),...
%             ClusterComparisonWrapper(@(a,b)(fuzzyComparison(a,b,'ari')), mVertPosMapAct),...
%             ClusterComparisonWrapper(@(a,b)(bmOverlapPosDist(a, mAdj, b, mAdj, 'nmi', 'discreteMemb', true)), mVertPosMapAct),...
%             ClusterComparisonWrapper(@(a,b)(bmOverlapPosDist(a, mAdj, b, mAdj, 'ami', 'discreteMemb', true)), mVertPosMapAct),...
%         };
%     else
%         cfValidMeasFunc = {};
%     end
        
 

    fImageDistanceFunc = ImagePurtity();        
%     fDistanceFunc = EucTriFact;  
       
    % initialise (this only tested for factorisation based algorithms)
    mInitMembership = fMemInitFunc.initMembership(mAdj, numPos);

    % initialise mImage
    mInitImage = fImageInitFunc.initImage(mAdj, mInitMembership);           

    algorRunNum = 1;
    
    

    if length(cAlgors) > 11
        cNewVarargin = cat(2, varargin, cAlgors{12:end});
    else
        cNewVarargin = varargin;
    end
        
    [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, runningTime, totalIterNum] =...
        runBMAlgor(mAdj, cAlgors{1}, cAlgors{2}, cAlgors{3}, cAlgors{4}, algorRunNum, cAlgors{6}, cAlgors{7}, cAlgors{8}, cAlgors{9}, cAlgors{10}, cAlgors{11},...
            fImageDistanceFunc, {}, 'initialMem', mInitMembership, 'initialImage', mInitImage, 'collectPerIterStats', true, cNewVarargin{:});        

    % output struct
    
    stResults.('bestObjective') = bestObjVal;
        
%     if sum(sum(~isreal(mBestMembership)))
%                     comparisonDistFuzzyAri = -1;
%         comparisonDistFuzzyJac = 0;
%         kktResidual = 10^6;
%         stResults.('fuzzyJac') = 0;
%         stResults.('kktResidual') = 10^6;
%     else
%                     comparisonDistFuzzyAri = fuzzyComparison(mVertPosMapAct, mBestMembership, 'ari');
%         comparisonDistFuzzyJac = fuzzyComparison(mVertPosMapAct, mBestMembership, 'jaccard');
%         kktResidual = fKKTDistanceFunc.residual(mAdj, mBestImage, mBestMembership);
%         stResults.('fuzzyJac') = fuzzyComparison(mVertPosMapAct, mBestMembership, 'jaccard');
%         stResults.('kktResidual') = fKKTDistanceFunc.residual(mAdj, mBestImage, mBestMembership);   
%     end
                
%     comparisonDistAmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'ami', 'discreteMemb', true);
%     comparisonDistNmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'nmi', 'discreteMemb', true);
%     stResults.('ami') = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'ami', 'discreteMemb', true);
%     stResults.('nmi') = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'nmi', 'discreteMemb', true);   
%                 comparisonDistAri = 1-bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'rand', 'discreteMemb', true);
                

                
    % seqiBloc encoding cost
    [~, bestAdjEncodingCost] = encodingCost(mAdj, mBestMembership);
    stResults.('bestAdjEncodingCost') = bestAdjEncodingCost;  
                
    % number of zro rows
    stResults.('zeroRowNum') = length(find(sum(mBestMembership,2) == 0));
                

                % store results
%                 vImageHist = histc(full(cmImage{c}(:)), 0:0.1:0.1*(imageBinNum-1));
%                 vPosHist = histc(full(cmMembership{c}(:)), 0:0.1:0.1*(posBinNum-1));

     
    figure;
    colormap(flipud(gray));
    imagesc(mBestMembership);
    
    figure;
    % plot image matrix
    colormap(flipud(gray));
    imagesc(mBestImage);
    
    mPosColour = [...
        0,0,1;...
        0,1,0;...
        1,0,0;...
        0.5, 0.5, 1;...
        0.5, 1, 0.5;...
        1, 0.5, 0.5;...
        ];
    
%     plotSingleBlockmodel(full(mAdj), discretise(mBestMembership), cAlgors{2}, false, false); 
    plotColouredSingleBlockmodel(full(mAdj), discretise(mBestMembership), mPosColour, cAlgors{2}, false, false); 
    
%     plotSingleBlockmodel(full(mBestMembership * mBestImage * mBestMembership'), discretise(mBestMembership), cAlgors{2}, false,false); 
    
    mBestMembership
    mBestImage    

                
    % display the various stats about each run of each algorithm 
    % (should have same number of positions whether single or hybrid, so
    % cccAlgors{a}{2}{3} = posNum).
%     if bCollectPerIterStats
%         plotInfo(sDataName, cAlgors{3}, mAdj, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmMembership, ccmImage);
%         mBestMembership
%         mBestImage
%     end

end % end of function



function plotInfo(sDataName, posNum, mAdj, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmMembership, ccmImage)
%
% Plots various stats about each algorithm results over the iterations of their
% runs.  One plot per run.
%

        
    % colour of lines
    cColour = {'b', 'r', 'y', 'g', 'm', 'c', 'k', 'b', 'r', 'y', 'g', 'm', 'c', 'k', 'b', 'r', 'y', 'g', 'm', 'c', 'k'};
    % plot the various things
    rowNum = 3;
    colNum = 5;
    % loop through each run
    for r = 1 : length(cvObjVal)
            figure;
            title(sDataName);
            iterNum = length(cvObjVal{r});
            subplot(rowNum, colNum, 1);
            plot(1:iterNum, cvObjVal{r});
            title('obj val');
            
            subplot(rowNum, colNum, 2);
            plot(1:iterNum, cvKKTResidual{r});
            title('KKT residual');
            
            subplot(rowNum, colNum, 3);
            plot(1:iterNum, cvImageDis{r});
            title('image residual');            
            
            subplot(rowNum, colNum, 4);
            plotyy(1:iterNum, ccvGroundComparison{r}{1}, 1:iterNum, ccvGroundComparison{r}{2});
            title('(fuzzy ground comparison');
            
            subplot(rowNum, colNum, 5);
            plotyy(1:iterNum, ccvGroundComparison{r}{3}, 1:iterNum, ccvGroundComparison{r}{4});
            title('(harden ground comparison');            
            
            % compute membership matrix stats
            mMemNNZNum = zeros(posNum, iterNum);
            mMemSum = zeros(posNum, iterNum);
            mMemAvg = zeros(posNum, iterNum);
            for c = 1 : iterNum
                for p = 1 : posNum
                    mMemNNZNum(p, c) = nnz(ccmMembership{r}{c}(:,p));
                    mMemSum(p, c) = sum(ccmMembership{r}{c}(:,p));
                    mMemAvg(p, c) = mean(ccmMembership{r}{c}(:,p));
                end
            end
            subplot(rowNum, colNum, 6);
            hold on;
            for p = 1 : posNum
                plot(1:iterNum, mMemNNZNum(p,:), cColour{p});
            end
            title('nnz of mMembership');
%             ylim([0, size(mAdj,1)]);
            hold off;
            
            subplot(rowNum, colNum, 7);
            hold on;
            for p = 1 : posNum
                plot(1:iterNum, mMemSum(p,:), cColour{p});
            end
            title('sum of mMembership');            
            hold off;         
            
            
            subplot(rowNum, colNum, 8);
            hold on;
            for p = 1 : posNum
                plot(1:iterNum, mMemAvg(p,:), cColour{p});
            end
            title('average of mMembership');            
            hold off;               
            
            
            % compute image matrix stats
            vImageNNZNum = zeros(1, iterNum);
            vImageSum = zeros(1, iterNum);
            for c = 1 : iterNum
                vImageNNZNum(c) = nnz(ccmImage{r}{c});
                vImageSum(c) = sum(sum(ccmImage{r}{c}));
            end
            subplot(rowNum, colNum, 9);
            plotyy(1:iterNum, vImageNNZNum, 1:iterNum, vImageSum);
            title('nnz and sum of mImage');
            
%             plotSingleBlockmodel(mAdj, mBestMembership, cccAlgors{a}{2}, true);
            % plot membership
            subplot(rowNum, colNum, 10);
            colormap(flipud(gray));
            imagesc(ccmMembership{r}{c});
            
            % plot image matrix
            subplot(rowNum, colNum, 11);
            colormap(flipud(gray));
            imagesc(ccmImage{r}{c});
            
            % plot approximation
            subplot(rowNum, colNum, 12);
            colormap(flipud(gray));
            imagesc(ccmMembership{r}{c} * ccmImage{r}{c} * ccmMembership{r}{c}');
            
            % plot adj matrix
        	subplot(rowNum, colNum, 13);
            colormap(flipud(gray));
            imagesc(mAdj);                
    end
end % end of function
