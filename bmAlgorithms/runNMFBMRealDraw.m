function runNMFBMRealDraw(sInBaseDir, sMatSuffix, sPartSuffix, sImageSuffix, numRun, vNumPositions, cDensityStats, sOutBaseDir, sOutPrefix,...
    bPosHeaderAct, bListFormatAct, vNum, bMatlabFormat, bSparse, bAddOneToIndex, bAddOneToPosVerts, bGraphHeader)
    
%
% Scans the input directories then runs the testing
%
% sInBaseDir - Base direction of input files, which includes the graph files,
% the known positions and the known image matrices.
% sMatSuffix - the file wildcard for graphs.
% sPartSuffix - the file wildcard for positions.
% sImageSuffix - the file wildcard for images.
% numRum - the number of runs to experiment over ([0,numRum]).  The appropriate
% input run files must exist.
% vNumPosition - vector of the number of generated positions test files to
% experiment over.
% cDensityStats - (proportion of dense blocks, density of dense blocks, std. of
% dense blocks, density of sparse blocks, std. of sparse blocks) to experiment
% over.
% sOutBaseDir - Base directory of output files.
% sOutPrefix - prefix of output files.
%
% bAddOne - whether to add 1 to vertex indices
%


    % self specify
    convEpsilon = 0.01;
    algorRunNum = 1;


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
    ccAlgors = { {'reichardtEqual', 'sa', 1, convEpsilon, algorRunNum, 'reichDeg'}...
            };              

% ccAlgors = { {'irm', 'gibbs', 1, convEpsilon, algorRunNum, 'irm'}...
%             };              
        
%    ccAlgors = {  {'matApprox', 'hardCProjGradDescM', 1, convEpsilon, algorRunNum, 'bmEuclidean'} };
        



    % get the input files
    % get the list of snapshot/matrices filenames
    sFileWildcard = sprintf('*.%s', sMatSuffix);
    tsMatFilenames = dir(fullfile(sInBaseDir, sFileWildcard));    
    
    % should only be one file
    assert(size(tsMatFilenames,1) == 1);

    sMatFile = tsMatFilenames(1).name;
    % get the associated position and image files
    sPosFile = strrep(sMatFile, sMatSuffix, sPartSuffix);
    sImageFile = strrep(sMatFile, sMatSuffix, sImageSuffix);

    % no graph or comparison out file
    [vBestObjVal, vBestMatApproxVal, vIdealEucObjVal, vIdealManObjVal, vRunningTime, vTotalIterTime] = ...
        runBMAlgorAndStore(fullfile(sInBaseDir, sMatFile), ...
            vNum, bSparse, bAddOneToIndex, bAddOneToPosVerts, bGraphHeader,...
                '', '', ccAlgors);
                
    assert(size(vBestObjVal,1) == size(ccAlgors,2));
                
    % store results
    for c = 1 : size(ccAlgors,2)
        cmResults{c}(r,:) = [vBestObjVal(c), vBestMatApproxVal(c), vIdealEucObjVal(c), vIdealManObjVal(c), vRunningTime(c), vTotalIterTime(c)];
    end     

        

    

    
    
end % end of function


   