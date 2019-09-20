function [vBestObjVal, vComparisonDistNmi, vComparisonDistAmi, vComparisonDistAri, vRunningTime, vTotalIterTime, vDegenerate] =...
    runBMAlgorAndCompare(sGraphFile, sActPosFile,...
        bPosHeaderAct, bListFormatAct, vNum, bMatlabFormat, bSparse, bAddOneToIndex, bAddOneToPosVerts, bGraphHeader,...
            sGraphFileOutPrefix, sComparisonOutFile, ccAlgors, varargin)
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
% vargin parameters relate to sigmoid parameters atm.
%
% sComparisonOutFile - name of the comparison output file.  Leave as '' if
% output not desired.
%


    % load the graph

    [mAdj] = loadMatlabGraph(sGraphFile, vNum, vNum, bMatlabFormat, bSparse, bAddOneToIndex, bGraphHeader);


    % load the ground truth positions
    [mVertPosMapAct, ~] = loadPositions(sActPosFile, bPosHeaderAct, bListFormatAct, vNum, bAddOneToPosVerts);


    % if empty, assign standard for ccAlgors
    if size(ccAlgors,2) == 0
        error('Need to specify ccAlgors');
    end

    % open approprite file for writing (APPEND mode)
    if size(sComparisonOutFile, 2) > 0
        fCompareResults = fopen(sComparisonOutFile,'a');
    end

    
    % comparison measures between ground truth and obtained positions
    cfValidMeasFunc = {...
%         ClusterComparisonWrapper(@(a,b)(bmOverlapPosDist(a, mAdj, b, mAdj, 'nmi', 'discreteMemb', true)), mVertPosMapAct),...
%         ClusterComparisonWrapper(@(a,b)(bmOverlapPosDist(a, mAdj, b, mAdj, 'ami', 'discreteMemb', true)), mVertPosMapAct),...
%         ClusterComparisonWrapper(@(a,b)(bmOverlapPosDist(a, mAdj, b, mAdj, 'rand', 'discreteMemb', true)), mVertPosMapAct),...
    };    
    
    
    vBestObjVal = zeros(size(ccAlgors,2), 1);
    vComparisonDistNmi = zeros(size(ccAlgors,2), 1);
    vComparisonDistAmi = zeros(size(ccAlgors,2), 1);
    vComparisonDistAri = zeros(size(ccAlgors,2), 1);
    vRunningTime = zeros(size(ccAlgors,2), 1);
    vTotalIterTime = zeros(size(ccAlgors,2), 1);
    vDegenerate = zeros(size(ccAlgors,2), 1);


    fImageDistanceFunc = EucImageDistance(zeros(size(mVertPosMapAct,2), size(mVertPosMapAct,2)));

    % run the various algorithms
    for a = 1 : size(ccAlgors,2)
        display(ccAlgors{a}{2});
        
        if length(ccAlgors{a}) > 11
            cNewVarargin = cat(2, varargin, ccAlgors{a}{12:end});
        else
            cNewVarargin = varargin;
        end

        

        [mBestImage, mBestMembership, bestObjVal, ~, ~, ~, ~, ~, ~, runningTime, totalIterNum] =...
            runBMAlgor(mAdj, ccAlgors{a}{1:11}, fImageDistanceFunc, cfValidMeasFunc, cNewVarargin{:});
        
    %     assert(size(mBestMembership,2) == posNum);

        % check if the number of solutions obtained is less than expected.  If so,
        % we ignore that result and log that the number of clusters was less than
        % actipicated
        vColSum = sum(mBestMembership, 1);
        if sum(vColSum == 0) > 0
            warning('Algorithm %s has produced a degenerate solution with less than the number of specified clusters', ccAlgors{a}{2});
            vDegenerate(a) = 1;
        end


%         % compare the positions     
%         comparisonDistVi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'vi');
%         comparisonDistAmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'ami');
%         % we use the version of nmi which is slower but does not require the
%         % two clusters to have the same number of positions.
%         % This is possible, e.g., BNMTF and some of the soft clustering approaches.
%         comparisonDistNmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'nmi');


        vBestObjVal(a) = bestObjVal;
        vComparisonDistNmi(a) = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'nmi', 'discreteMemb', true);
        vComparisonDistAmi(a) = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'ami', 'discreteMemb', true);
        vComparisonDistAri(a) = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'rand', 'discreteMemb', true);
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
            fprintf(fCompareResults, '%s, %s, %d, %6f, %6f, %9f, %9f\n',...
                ccAlgors{a}{1}, ccAlgors{a}{2}, ccAlgors{a}{3}, bestObjVal, idealObjVal, comparisonDist, runningTime);
        end % end of writing out to file
    end % end of outer for
        
    if size(sComparisonOutFile, 2) > 0
        fclose(fCompareResults);
    end


end % end of function




