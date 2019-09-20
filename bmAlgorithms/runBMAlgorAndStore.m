function [vBestObjVal, vBestMatApproxVal, vIdealEucObjVal, vIdealManObjVal, vRunningTime, vTotalIterTime] = runBMAlgorAndStore(sGraphFile, ...
    vNum, bMatlabFormat, bSparse, bAddOneToIndex, bGraphHeader, sGraphFileOutPrefix, sComparisonOutFile, ccAlgors, bDrawBlockmodel)
%
% Loads a ground truth BM, then runs the graph on various algorithms, then
% compares the found positions against the ground truth, then appends the
% results into the results file.
%
% bPosHeaderAct, bListFormatAct, bAddOneToPosVerts related to the loading
% of the ground truth for positions.
% bMatlabFormat, bSparse, bAddOneToIndex, bGraphHeader related to graph
% loading.
%
% sComparisonOutFile - name of the comparison output file.  Leave as '' if
% output not desired.
%



    % self specify (default)
    convEpsilon = 0.01;
    algorRunNum = 1;
    runNum = 1;


    % algorithm specifications (set default of one position, until they are set
%     % within the loops
%     ccAlgors = { {'matApprox', 'hardChardM', 1, convEpsilon, algorRunNum, 'euclidean'},...
%             {'matApprox', 'hardCkmeansM', 1, convEpsilon, algorRunNum, 'euclidean'},...
%             {'matApprox', 'hardCbolongM', 1, convEpsilon, algorRunNum, 'euclidean'},...
%             {'matApprox', 'hardCProjGradDescM', 1, convEpsilon, algorRunNum, 'euclidean'},...
%             {'matApprox', 'hardCCoordDescM', 1, convEpsilon, algorRunNum, 'euclidean'}...
%             {'matApprox', 'hardCCoordDescExactM', 1, convEpsilon, algorRunNum, 'euclidean'}...
%             {'matApprox', 'hardCProjGradDescM', 1, convEpsilon, algorRunNum, 'bmEuclidean'},...
%             {'matApprox', 'hardCCoordDescM', 1, convEpsilon, algorRunNum, 'bmEuclidean'}...
%             };
        
%     ccAlgors = { {'hardSa', 'sa', 1, convEpsilon, algorRunNum, 'bmEuclidean'}...
%             };       
%     ccAlgors = { {'reichardtEqual', 'sa', 1, convEpsilon, algorRunNum, 'reichDeg'}...
%             };              

% ccAlgors = { {'irm', 'gibbs', 1, convEpsilon, algorRunNum, 'irm'}...
%             };              
        
%    ccAlgors = {  {'matApprox', 'hardCProjGradDescM', 1, convEpsilon, algorRunNum, 'bmEuclidean'} };


% 
%     % get the input files
%     % get the list of snapshot/matrices filenames
%     sFileWildcard = sprintf('*.%s', sMatSuffix);
%     tsMatFilenames = dir(fullfile(sInBaseDir, sFileWildcard));    
%     
%     % should only be one file
%     assert(size(tsMatFilenames,1) == 1);
% 
%     sMatFile = tsMatFilenames(1).name;
% 
%     % no graph or comparison out file
%     [vBestObjVal, vBestMatApproxVal, vIdealEucObjVal, vIdealManObjVal, vRunningTime, vTotalIterTime] = ...
%         runBMAlgorAndStore(fullfile(sInBaseDir, sMatFile), ...
%             vNum, bSparse, bAddOneToIndex, bAddOneToPosVerts, bGraphHeader,...
%                 '', '', ccAlgors);
%                 
%     assert(size(vBestObjVal,1) == size(ccAlgors,2));
%                 
%     % store results
%     for c = 1 : size(ccAlgors,2)
%         cmResults{c}(r,:) = [vBestObjVal(c), vBestMatApproxVal(c), vIdealEucObjVal(c), vIdealManObjVal(c), vRunningTime(c), vTotalIterTime(c)];
%     end     






    % load the graph
    [mAdj] = loadMatlabGraph(sGraphFile, vNum, vNum, bMatlabFormat, bSparse, bAddOneToIndex, bGraphHeader);


    % open approprite file for writing (APPEND mode)
    if size(sComparisonOutFile, 2) > 0
        fCompareResults = fopen(sComparisonOutFile,'a');
    end

    vBestObjVal = zeros(size(ccAlgors,2), 1);
    vBestMatApproxVal = zeros(size(ccAlgors,2), 1);
    vIdealEucObjVal = zeros(size(ccAlgors,2), 1);
    vIdealManObjVal = zeros(size(ccAlgors,2), 1);
    vRunningTime = zeros(size(ccAlgors,2), 1);
    vTotalIterTime = zeros(size(ccAlgors,2), 1);


    % run the various algorithms
    for a = 1 : size(ccAlgors,2)
        display(ccAlgors{a}{2});
    
        [mBestImage, mBestMembership, bestObjVal, bestMatApproxVal, idealEucObjVal, idealManObjVal, runningTime, totalIterNum] = runBMAlgor(mAdj,...
            ccAlgors{a}{1}, ccAlgors{a}{2}, ccAlgors{a}{3}, ccAlgors{a}{4}, ccAlgors{a}{5}, ccAlgors{a}{6}, ccAlgors{a}{7});

    
        vBestObjVal(a) = bestObjVal;
        vBestMatApproxVal(a) = bestMatApproxVal;
        vIdealEucObjVal(a) = idealEucObjVal;
        vIdealManObjVal(a) = idealManObjVal;
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
    
        if size(sComparisonOutFile, 2) > 0
            % write out the comparison results
            fprintf(fCompareResults, '%s, %s, %d, %6f, %6f, %6f %6f\n',...
                ccAlgors{a}{1}, ccAlgors{a}{2}, ccAlgors{a}{3}, bestObjVal, idealEucObjVal, idealManObjVal, runningTime);
        end % end of writing out to file
        
        
        % see if we should draw the blockmodel
        if bDrawBlockmodel
           plotSingleBlockmodel(mAdj, mBestMembership, ccAlgors{a}{2}); 
        end
    end
        
    if size(sComparisonOutFile, 2) > 0
        fclose(fCompareResults);
    end







end % end of function




