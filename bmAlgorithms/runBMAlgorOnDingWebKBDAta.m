function runBMAlgorOnDingWebKBDAta(cData,...
        numRun, sOutBaseDir, sOutPrefix, varargin)
            
    
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
    addParameter(inParser, 'drawBlockmodel', false);
    addParameter(inParser, 'convergEsp', 0.01);
    addParameter(inParser, 'algorRunNum', 1);

    parse(inParser, varargin{:});
    
    bDiscretiseMembership = inParser.Results.discreteMemb;
    bColNorm = inParser.Results.normaliseCol;
    bDrawBlockmodel = inParser.Results.drawBlockmodel;
    convEpsilon = inParser.Results.convergEsp;
    algorRunNum = inParser.Results.algorRunNum;
    
    % get remaining arguments
    cFields = fieldnames(inParser.Unmatched);
    cRemainArgs = cell(1, 2*numel(cFields));
    for i = 1 : numel(cFields)
        cRemainArgs{i} = cFields(i);
        cRemainArgs{i+1} = inParser.Unmatched.(cFields{i});
    end

    
    % self specify
    initPosNum = 1; %this will be updated with the actual posNum when the ground truths are loaded

    bDonotDiscretise = false;
    

    % algorithm specifications (set default of one position, until they are set
    
    % missing bnmtf
    %         {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', bDiscretiseMembership},...
    ccAlgors = {...
%         {'matApprox', 'hardM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
%         {'matApprox', 'kmeansM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
%         {'matApprox', 'bolongM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
%         {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCIncr', bDonotDiscretise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCIncr', bDonotDiscretise},...
%         {'reichardt', 'saHamilM', 1, convEpsilon, algorRunNum, 'reichEqual', 'saHamilM', bDonotDiscretise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise},...
%         {'reichardt', 'saHamilM', 1, convEpsilon, algorRunNum, 'reichDeg', 'saHamilM', bDonotDiscretise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'hardCIncrAdjDeg', bDonotDiscretise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'hardCIncrAdjDeg', bDonotDiscretise},...
%         {'hardSa', 'saHardM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'saHardC', bDonotDiscretise},...
        {'matFact', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', bDiscretiseMembership, bColNorm},...
        {'matFact', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', bDiscretiseMembership, bColNorm},...
%             {'irm', 'irm', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'irm', bDiscretiseMembership, bColNorm},...
%             {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', bDiscretiseMembership, bColNorm},...
%             {'bkn', 'bkn', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'bkn', bDiscretiseMembership, bColNorm},...        
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCGradDescAdjEqual', bDiscretiseMembership, bColNorm},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCCoordDescAdjEqual', bDiscretiseMembership, bColNorm},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCGradDescAdjDeg', bDiscretiseMembership, bColNorm},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCCoordDescAdjDeg', bDiscretiseMembership, bColNorm},...
        };

% 
%     cData = cell(4,3);
%     cData{1,1} = double(mAdjCornell);
%     cData{2,1} = double(mAdjTexas);
%     cData{3,1} = double(mAdjWashington);
%     cData{4,1} = double(mAdjWisconsin);
%     cData{1,2} = mPosActCornell;
%     cData{2,2} = mPosActTexas;
%     cData{3,2} = mPosActWashington;
%     cData{4,2} = mPosActWisconsin;
%     cData{1,3} = 'cornell';
%     cData{2,3} = 'texas';
%     cData{3,3} = 'washington';
%     cData{4,3} = 'wisconsin';
    
    
    
    cfOutResults = cell(size(cData,1),1);
        
    % open a result file per dataset
    for c = 1 : size(cData,2)
        sOutResultFile = sprintf('%s_%s%s', sOutPrefix, cData{c}{3},...
            '.result.csv');  
        cfOutResults{c} = fopen(fullfile(sOutBaseDir, sOutResultFile), 'a');
    end
            
    performRuns(cfOutResults, ccAlgors, cData, numRun, bDrawBlockmodel);
    
    % close the files
    for c = 1 : size(cData,2)
        fclose(cfOutResults{c});
    end
    
        
end % end of function



function performRuns(cfOutResults, ccAlgors, cData, numRun, bDrawBlockmodel)

%     numPos = 7;
    
    % loop thru the datasets
    for d = 1 : size(cData,2)
        mAdj = cData{d}{1};
        mVertPosMapAct = cData{d}{2};
        sDataName = cData{d}{3};
        numPos = length(unique(mVertPosMapAct(:,2)));
        
        display(sDataName);

        % number of results per run per algor
        resultNum = 9;
        cmResults = cell(size(ccAlgors,2),1);
        for c = 1 : size(ccAlgors,2)
            cmResults{c} = zeros(numRun, resultNum);
            ccAlgors{c}{3} = numPos;
        end
            

        
            % loop through the runs
            for r = 1 : numRun          
                % compare
                [vBestObjVal, vBestMatApproxVal, vIdealEucObjVal, vIdealManObjVal, vComparisonDistVi, vComparisonDistAmi, vComparisonDistNmi, vRunningTime, vTotalIterTime] = ...
                    runAndCompare(mAdj, mVertPosMapAct, ccAlgors, bDrawBlockmodel, '');
                
                assert(size(vBestObjVal,1) == size(ccAlgors,2));
    
                % store results
                for c = 1 : size(ccAlgors,2)
                    cmResults{c}(r,:) = [vBestObjVal(c), vBestMatApproxVal(c), vIdealEucObjVal(c), vIdealManObjVal(c), vComparisonDistVi(c), vComparisonDistAmi(c), vComparisonDistNmi(c), vRunningTime(c), vTotalIterTime(c)];
                end

            end % end of loop through numRum

            
        % compute average
        for c = 1 : size(ccAlgors,2)
            vAvg = mean(cmResults{c},1);
            if numRun > 1
                vStd = std(cmResults{c},1);
            else
                vStd = zeros(1, size(cmResults{c}, 2));
            end
            sAlgorname = sprintf('%s-%s-%s-%s', ccAlgors{c}{1}, ccAlgors{c}{2}, ccAlgors{c}{6}, ccAlgors{c}{7});
%             disp(sprintf( '%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n', vAvg, vStd));
            fprintf(cfOutResults{d}, '%s, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n', sAlgorname, vAvg, vStd);  
        end
    
    end % end of loop thru datasets


end % end of function


function [vBestObjVal, vBestMatApproxVal, vIdealEucObjVal, vIdealManObjVal, vComparisonDistVi, vComparisonDistAmi, vComparisonDistNmi, vRunningTime, vTotalIterTime] =...
    runAndCompare(mAdj, mVertPosMapAct, ccAlgors, bDrawBlockmodel, sGraphFileOutPrefix)




    vBestObjVal = zeros(size(ccAlgors,2), 1);
    vBestMatApproxVal = zeros(size(ccAlgors,2), 1);
    vIdealEucObjVal = zeros(size(ccAlgors,2), 1);
    vIdealManObjVal = zeros(size(ccAlgors,2), 1);
    vComparisonDistVi = zeros(size(ccAlgors,2), 1);
    vComparisonDistAmi = zeros(size(ccAlgors,2), 1);
    vComparisonDistNmi = zeros(size(ccAlgors,2), 1);
    vRunningTime = zeros(size(ccAlgors,2), 1);
    vTotalIterTime = zeros(size(ccAlgors,2), 1);



    % run the various algorithms
    for a = 1 : size(ccAlgors,2)
        display(ccAlgors{a}{2});

        
        [mBestImage, mBestMembership, bestObjVal, bestMatApproxVal, idealEucObjVal, idealManObjVal, runningTime, totalIterNum] = runBMAlgor(mAdj,...
            ccAlgors{a}{1}, ccAlgors{a}{2}, ccAlgors{a}{3}, ccAlgors{a}{4}, ccAlgors{a}{5}, ccAlgors{a}{6}, ccAlgors{a}{7}, ccAlgors{a}{8}, ccAlgors{a}{9});
    
        comparisonDistVi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'vi');
        comparisonDistAmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'ami');
        comparisonDistNmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'nmi');

    
        vBestObjVal(a) = bestObjVal;
        vBestMatApproxVal(a) = bestMatApproxVal;
        vIdealEucObjVal(a) = idealEucObjVal;
        vIdealManObjVal(a) = idealManObjVal;
        vComparisonDistVi(a) = comparisonDistVi;
        vComparisonDistAmi(a) = comparisonDistAmi;
        vComparisonDistNmi(a) = comparisonDistNmi;
        vRunningTime(a) = runningTime;
        vTotalIterTime(a) = totalIterNum;
    
    
        if size(sGraphFileOutPrefix, 2) > 0
            sMemOutFile = sprintf('%s_%s_%s_%d%s', sGraphFileOutPrefix, ccAlgors{a}{1}, ccAlgors{a}{2}, ccAlgors{a}{3}, '.part.csv');
    
            % write out mBestMembership
            fMemOut = fopen(sMemOutFile, 'w');
            % find the members in each position and output to file
            for p = 1 : size(mBestMembership, 2)
                vNZ = find(mBestMembership(:,p));
                vNZ = vNZ - 1;
                % we need ot subtract 1 to get it back to 0 based indexing
                fprintf(fMemOut, '%d,',vNZ(1:end-1));
                fprintf(fMemOut, '%d\n', vNZ(end));
            end
            fclose(fMemOut);
    
            % write out mBestImage
            if issparse(mBestImage)
                mBestImage = full(mBestImage);
            end
    
            sImageOutfile = sprintf('%s_%s_%s_%d%s', sGraphFileOutPrefix, ccAlgors{a}{1}, ccAlgors{a}{2}, ccAlgors{a}{3}, '.image.csv');
    
            csvwrite(sImageOutfile, mBestImage);
        end
    

               
        % see if we should draw the blockmodel
        if bDrawBlockmodel
           plotSingleBlockmodel(mAdj, mBestMembership, ccAlgors{a}{2}); 
        end
    end
        
%     if size(sComparisonOutFile, 2) > 0
%         fclose(fCompareResults);
%     end







end % end of function




