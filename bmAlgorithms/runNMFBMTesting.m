function runNMFBMTesting(sInBaseDir, sMatSuffix, sPartSuffix, sImageSuffix, datasetNumPerStat, runsPerDataset, vNumPositions,...
    cDensityStats, sOutBaseDir, sOutPrefix, varargin)
    
%
% Used for testing the effect of synthetic datasets (density, noise) on
% factorisation algorithms.
%
% Scans the input directories then runs the testing.  For image matrix with
% sigmoid loss term.
%
% sInBaseDir - Base direction of input files, which includes the graph files,
% the known positions and the known image matrices.
% sMatSuffix - the file wildcard for graphs.
% sPartSuffix - the file wildcard for positions.
% datasetNumPerStat - the number of runs to experiment over ([0,numRum]).  The appropriate
% input run files must exist.
% runsPerDataset - the number of runs to do per dataset.
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
% vSigmoidWeight, vSigmoidGamma - extra parameters to do with the bmEucildean
% distance function and variants.  
%
% @author: Jeffrey Chan, 2013
%


    % parse arguments
    inParser = inputParser;
    inParser.KeepUnmatched = true;

    % add parameters and set default valuse
    addParameter(inParser, 'discreteMemb', false);
    addParameter(inParser, 'normaliseCol', false);
    addParameter(inParser, 'convergeEsp', 0.0001);
    addParameter(inParser, 'vSigmoidWeight', [100000]);
    addParameter(inParser, 'vHardRelaxWeight', [1]);
    addParameter(inParser, 'algorRunNum', 1);

    parse(inParser, varargin{:});
    
    bDiscretiseMembership = inParser.Results.discreteMemb;
    bColNormalise = inParser.Results.normaliseCol;
    convEpsilon = inParser.Results.convergeEsp;
    vSigmoidWeight = inParser.Results.vSigmoidWeight;
    vHardRelaxWeight = inParser.Results.vHardRelaxWeight;
    algorRunNum = inParser.Results.algorRunNum;
    
    initPosNum = 1; %this will be updated with the actual posNum when the ground truths are loaded
    
    % for algorithms that should not be discretised (because they are hard
    % algoritms.
    bDonotDiscretise = false;
    bDonotColNormalise = false;
    % for others, it depends on bDiscretiseMembership
    


%    
%     ccAlgors = {  
%         {'matApprox', 'hardM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...

    
%         {'matApprox', 'hardM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'bolongM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCIncr', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCIncr', bDonotDiscretise, bDonotColNormalise},...        
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjDeg', 'hardCIncrAdjDeg', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjDeg', 'hardCIncrAdjDeg', bDonotDiscretise, bDonotColNormalise},...        
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'hardCIncrAdjDeg', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'hardCIncrAdjDeg', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'bolongM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCGradDescAdjEqual', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCCoordDescAdjEqual', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCGradDescAdjDeg', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCCoordDescAdjDeg', bDiscretiseMembership, bColNormalise},...        
%         {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', bDiscretiseMembership, bDonotColNormalise, bColNormalise},...
%         };



    ccAlgors = {
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bDonotColNormalise},...
        {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bDonotColNormalise},...
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bDonotColNormalise},...
        {'matApprox', 'coordDescBMM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bDonotColNormalise},...
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'softCMultHardRowConstraintAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bDonotColNormalise},...
        {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'softCMultHardRowConstraintAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bDonotColNormalise},...
        {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCMultHardRowConstraintAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bDonotColNormalise},...
        {'matApprox', 'coordDescBMM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCMultHardRowConstraintAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bDonotColNormalise},...
%         {'matApprox', 'bolongM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', 'randomInit', 'randomInit', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', 'randomInit', 'randomInit', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', 'randomInit', 'randomInit', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', 'randomInit', 'randomInit', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCIncr', 'randomInit', 'randomInit', bDonotDiscretise, bDonotColNormalise, 'sigmoidWeight', 0.1},...
%         {'matApprox', 'coordDescBMM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCIncr', 'randomInit', 'randomInit', bDonotDiscretise, bDonotColNormalise, 'sigmoidWeight', 0.1},...        
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'hardCIncrAdjEqual', 'randomInit', 'randomInit', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'hardCIncrAdjEqual', 'randomInit', 'randomInit', bDonotDiscretise, bDonotColNormalise},...
%         {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', true, bColNormalise, 'runIterationNum', 400, 'rowNormaliseAdhoc', false},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', true, bColNormalise, 'runIterationNum', 400, 'rowNormaliseAdhoc', false},...        
%         {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', 'randomInit', 'randomInit', true, bColNormalise, 'runIterationNum', 100},...        
%         {'matApprox', 'bolongM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', bDiscretiseMembership},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', bDiscretiseMembership},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCGradDescAdjEqual', bDiscretiseMembership},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCCoordDescAdjEqual', bDiscretiseMembership},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCGradDescAdjDeg', bDiscretiseMembership}        {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCCoordDescAdjDeg', bDiscretiseMembership},...        
%         {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', bDiscretiseMembership},...
%         {'matApprox', 'bolongM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCGradDescAdjEqual', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCCoordDescAdjEqual', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCGradDescAdjDeg', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCCoordDescAdjDeg', bDiscretiseMembership, bColNormalise},...  
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCGradDescAdjEqual', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCCoordDescAdjEqual', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCGradDescAdjDeg', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjDeg', 'softCCoordDescAdjDeg', bDiscretiseMembership, bColNormalise},...         
%         {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', bDiscretiseMembership, bColNormalise}
        };

    
% 
% 
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'softCGradDescAdjEqual', bDiscretiseMembership},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCGradDesc', bDiscretiseMembership},...
%         {'matApprox', 'coordDescBMM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCGradDescAdjEqual', bDiscretiseMembership},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'softCCoordDescAdjEqual', bDiscretiseMembership},...  
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCCoordDesc', bDiscretiseMembership},...  
%         {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', bDonotDiscretise},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCIncr', bDonotDiscretise},...        
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise},...        
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'hardCIncrAdjEqual', bDonotDiscretise},...


    
%     ccAlgors = {...
% % %         {'mmsbm', 'mmsbmM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCmmsbm', 'randomInit', 'randomInit', bDiscretiseMembership},...                
% {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon+1, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'runIterationNum', 400, 'rowNormaliseAdhoc', true},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon+1, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'runIterationNum', 400, 'rowNormaliseAdhoc', true},...
%         {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'runIterationNum', 400, 'rowNormaliseAdhoc', false},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'runIterationNum', 400, 'rowNormaliseAdhoc', false},...        
%         {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'runIterationNum', 100},...        
%         {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'sigmoidWeight', 0.0, 'runIterationNum', 100},...
%         {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'sigmoidWeight', 0.0, 'runIterationNum', 100},...
%         {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraintRepara', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'sigmoidWeight', 0.0, 'runIterationNum', 100},...
%         {'matApprox', 'multLangM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'sigmoidWeight', 0.0, 'runIterationNum', 100},...
%         {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'sigmoidWeight', 0.0, 'runIterationNum', 100, 'regWeight', 10},...
%         {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'runIterationNum', 100, 'regWeight', 10},...        
%         {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'runIterationNum', 100},... 
%         {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'runIterationNum', 100},...
%         {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultHardRowConstraintRepara', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise, 'runIterationNum', 100},...              
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCGradDescAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise},...
%         {'matApprox', 'coordDescBMM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'softCCoordDescAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise},...           
%         {'matApprox', 'coordDescBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclideanAdjEqual', 'softCCoordDescAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise},...           
%         {'bkn', 'bknM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCbkn', 'randomInit', 'randomInit', bDiscretiseMembership, bColNormalise},...        
%         };
    





    


    % dataset file info
    %cDatasetInfo = {bPosHeaderAct, bListFormatAct, bMatlabFormat, bSparse, bAddOneToIndex, bAddOneToPosVerts, bGraphHeader};
    cDatasetInfo = {false, true, true, false, false, true, false};


    % if no density parameters specified, we just run the algorithm
    if isempty(cDensityStats)
        cfOutResults = cell(size(ccAlgors,2),1);
        
        % open a result file per algorithm
        for c = 1 : size(ccAlgors,2)
            sOutResultFile = sprintf('%s_%s_%s_%s_%.2f_%d%s', sOutPrefix, ccAlgors{c}{2}, ccAlgors{c}{7}, ccAlgors{c}{6}, ccAlgors{c}{4}, ccAlgors{c}{5},...
                '.part.csv')  
            cfOutResults{c} = fopen(fullfile(sOutBaseDir, sOutResultFile), 'a');
        end
        
        loopThruPositions(cfOutResults, vNumPositions(1), vNumPositions(2), ccAlgors, 'sparse',...
            sInBaseDir, sMatSuffix, sPartSuffix, datasetNumPerStat, runsPerDataset,...
            cDatasetInfo{1}, cDatasetInfo{2}, vNum, cDatasetInfo{3}, cDatasetInfo{4}, cDatasetInfo{5}, cDatasetInfo{6}, cDatasetInfo{7})         
    end
    
    
    sImageInit = 'randomInit';
    sMemInit = 'randomInit';    
        
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
    

    
    % construct outfiles
    cfResultsFiles = cell(1, size(ccAlgors,2));
    for c = 1 : size(ccAlgors,2)
        sOutResultFile = sprintf('%s_%s_%s_%s_%.2f_%d%s', sOutPrefix, ccAlgors{c}{2}, ccAlgors{c}{7}, ccAlgors{c}{6}, ccAlgors{c}{4}, ccAlgors{c}{5},...
            '.results.csv');
        cfResultsFiles{c} = fopen(fullfile(sOutBaseDir, sOutResultFile), 'a'); 
    end
    
            
    % loop through the stats
    for s = 1 : length(cDensityStats)
        vStats = cDensityStats{s};
        fprintf('Doing stats %.2f, %.2f, %.2f, %.2f, %.2f, %d.\n', vStats);
            
        vNum = vStats(6);
        
        % loop through the known number of positions in the dataset
        for p = vNumPositions(1) : vNumPositions(2)
    %       fprintf('Doing position %d.\n', p);

            % number of results per run per (runPerDataset) per algor
            resultNum = 10;
            degenCol = 10;
            vRestCol = 1:9;
            nmiCol = 4;            

            % initialise cmResults
            cccmResults = cell(1, length(vSigmoidWeight));
            for sig = 1 : length(vSigmoidWeight)
                cccmResults{sig} = cell(1, length(vHardRelaxWeight));
                for rel = 1 : length(vHardRelaxWeight)
                    cccmResults{sig}{rel} = cell(1, size(ccAlgors,2));
                    for c = 1 : size(ccAlgors,2)
        %                 cmResults{c} = zeros(datasetNumPerStat * runsPerDataset, resultNum);
                        cccmResults{sig}{rel}{c} = zeros(datasetNumPerStat, resultNum);
                        ccAlgors{c}{3} = p;
                    end                         
                end
            end
        
            % loop through the different datasets per stat set
            for r = 1 : datasetNumPerStat
    %                             fprintf('run %d\n', r);
                % get the list of snapshot/matrices filenames
                sFileWildcard = sprintf('*%d-%d-%d*%.2f-%.2f-%.2f-%.2f-%.2f*%s', r, vNum, p, ...
                    vStats(1), vStats(2), vStats(3), vStats(4), vStats(5), sMatSuffix);
                fprintf('%s\n', sFileWildcard);
                tsMatFilenames = dir(fullfile(sInBaseDir, sFileWildcard));

                % should only be one file
                if size(tsMatFilenames,1) ~= 1
                    error('runNMFBMTesting:blah', 'There should only be one file returned from filewildcard.');
                end       

                sMatFile = tsMatFilenames(1).name;
                % get the associated position and image files
                sPosFile = strrep(sMatFile, sMatSuffix, sPartSuffix);       
                sImageFile = strrep(sMatFile, sMatSuffix, sImageSuffix);

                
                % load the graph
                [mAdj] = loadMatlabGraph(fullfile(sInBaseDir, sMatFile), vNum, vNum, cDatasetInfo{3}, cDatasetInfo{4}, cDatasetInfo{5}, cDatasetInfo{7});

                % load the ground truth positions
                [mVertPosMapAct, ~] = loadPositions(fullfile(sInBaseDir, sPosFile), cDatasetInfo{1}, cDatasetInfo{2}, vNum, cDatasetInfo{6});                
                
                [mImageAct] = loadImage(fullfile(sInBaseDir, sImageFile), size(mVertPosMapAct,2));
                
       
                
                
                % loop through sigmoid parameters
                for w = 1 : length(vSigmoidWeight)               
        %                     fprintf('Doing sigmoid weight %.2f.\n', vSigmoidWeight(w));
        
                    for phi = 1 : length(vHardRelaxWeight)
                        
                        cCurrDataResult = cell(1, size(ccAlgors,2));
                        for c = 1 : size(ccAlgors,2)
            %                 cmResults{c} = zeros(datasetNumPerStat * runsPerDataset, resultNum);
                            cCurrDataResult{c} = zeros(runsPerDataset, resultNum);
                        end                                    
                        
                        for i = 1 : runsPerDataset
    %                     fprintf('runsPerDataset %d\n', i);   

                            % initialise (this only tested for factorisation based algorithms)
                            mInitMembership = fMemInitFunc.initMembership(mAdj, size(mVertPosMapAct,2));

                            % initialise mImage
                            mInitImage = fImageInitFunc.initImage(mAdj, mInitMembership);               

                            for c = 1 : size(ccAlgors,2)
                                [bestObjVal,  bestImageDis, comparisonDistNmi, comparisonDistAmi, comparisonDistAri, comparisonDistFuzzyAri, comparisonDistFuzzyJac, runningTime, totalIterTime, degenSolution] = ...
                                    runAndCompare(mAdj, mVertPosMapAct, mImageAct, ccAlgors{c},...
                                        'sigmoidWeight', vSigmoidWeight(w), 'penaltyWeight', vHardRelaxWeight(phi), 'initialMem', mInitMembership, 'initialImage', mInitImage, varargin{:});

                                % store results
%                                 cmResults{c}( (r-1)*runsPerDataset + i, :) = [bestObjVal, bestImageDis, comparisonDistNmi, comparisonDistAmi, comparisonDistAri,  comparisonDistFuzzyAri, comparisonDistFuzzyJac, runningTime, totalIterTime, degenSolution];
                                cCurrDataResult{c}(i, :) = [bestObjVal, bestImageDis, comparisonDistAmi, comparisonDistNmi, comparisonDistAri,  comparisonDistFuzzyAri, comparisonDistFuzzyJac, runningTime, totalIterTime, degenSolution];
                            end % end of algorithms loop
                        end % end of loop through runsPerDataset
      
                        % find maximum

                        for c = 1 : size(ccAlgors,2) 
                            % find max
                            vRows = find(cCurrDataResult{c}(:,degenCol) == 0);
                            [~, maxIndex] = max(cCurrDataResult{c}(vRows,nmiCol));
                            cccmResults{w}{phi}{c}(r,:) = cCurrDataResult{c}(vRows(maxIndex),:);
                        end
                  
                        
                    end % end of loop through vHardRelaxWeight
                    
                end % end of loop through sigmoid
                
            end % end of looop through datasetNumPerStat
            

            % compute average
            for sig = 1 : length(vSigmoidWeight)    
                for rel = 1 : length(vHardRelaxWeight)
                    for c = 1 : size(ccAlgors,2)
        %                 vRows = find(cmResults{c}(:,degenCol) == 0);
                        vAvg = mean(cccmResults{sig}{rel}{c}(:,vRestCol),1);
                        if length(vRows) <= 1
                            vStd = zeros(1, length(vRestCol));
                        else
                            vStd = std(cccmResults{sig}{rel}{c}(:,vRestCol),1);
                        end
                        degenNum = sum(cccmResults{sig}{rel}{c}(:,degenCol) == 0);
                        fprintf('%d, %d, %.5f, %.5f, %.2f, %.2f, %.2f, %.2f, %.2f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %d\n', vNum, p, vSigmoidWeight(sig), vHardRelaxWeight(rel), vStats(1:5), vAvg, vStd, degenNum);                          
                        fprintf(cfResultsFiles{c}, '%d, %d, %.5f, %.5f, %.2f, %.2f, %.2f, %.2f, %.2f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %d\n', vNum, p, vSigmoidWeight(sig), vHardRelaxWeight(rel), vStats(1:5), vAvg, vStd, degenNum);                          
                    end                         
                end
            end            
            
               
                
        end % loop through number of positions
                
    end % loop through density stats
        
    for c = 1 : size(ccAlgors,2)
        fclose(cfResultsFiles{c});
    end
    
    % make sure all are closed
    fclose('all');
    
    
%     end % loop thru algorithms    
    
    
%     % open a result file per algorithm
%     for c = 1 : size(ccAlgors,2)
%         sOutResultFile = sprintf('%s_%s_%s_%s_%.2f_%d%s', sOutPrefix, ccAlgors{c}{2}, ccAlgors{c}{7}, ccAlgors{c}{6}, ccAlgors{c}{4}, ccAlgors{c}{5},...
%             '.results.csv')
%             
%         fOutResults = fopen(fullfile(sOutBaseDir, sOutResultFile), 'a');
%     
%         
%         % loop through the stats
%         for s = 1 : length(cDensityStats)
%             vStats = cDensityStats{s};
%             fprintf('Doing stats %.2f, %.2f, %.2f, %.2f, %.2f, %d.\n', vStats);
%             
%             vNum = vStats(6);
%             
%             % loop through the number of positions
%             for p = vNumPositions(1) : vNumPositions(2)
% %                 fprintf('Doing position %d.\n', p);
% 
%                 ccAlgors{c}{3} = p;
% 
% 
%                 % loop through sigmoid parameters
%                 for w = 1 : length(vSigmoidWeight)
% %                     fprintf('Doing sigmoid weight %.2f.\n', vSigmoidWeight(w));
%                     for g = 1 : length(vSigmoidGamma)
% %                         fprintf('Doing sigmoid gamma %.2f.\n', vSigmoidGamma(g));
%                         % number of results per run per (runPerDataset) per algor
%                         resultNum = 7;
%                         degenCol = 7;
%                         vRestCol = [1:6];
%                         mResults = zeros(numRun * runsPerDataset, resultNum);
%                         
%                         % loop through the runs
%                         for r = 1 : numRun          
% %                             fprintf('run %d\n', r);
%                             % get the list of snapshot/matrices filenames
%                             sFileWildcard = sprintf('*%d-%d-%d*%.2f-%.2f-%.2f-%.2f-%.2f*%s', r, vNum, p, ...
%                                 vStats(1), vStats(2), vStats(3), vStats(4), vStats(5), sMatSuffix);
%                             fprintf('%s\n', sFileWildcard);
%                             tsMatFilenames = dir(fullfile(sInBaseDir, sFileWildcard));
%                 
%                             % should only be one file
%                             if size(tsMatFilenames,1) ~= 1
%                                 error('runNMFBMTesting:blah', 'There should only be one file returned from filewildcard.');
%                             end
% 
% 
%                             sMatFile = tsMatFilenames(1).name;
%                             % get the associated position and image files
%                             sPosFile = strrep(sMatFile, sMatSuffix, sPartSuffix);
%                     
%                             for i = 1 : runsPerDataset
% %                                 fprintf('runsPerDataset %d\n', i);                    
%                                 [vBestObjVal, vComparisonDistNmi, vComparisonDistAmi, vComparisonDistAri, vRunningTime, vTotalIterTime, vDegenSolution] = ...
%                                     runBMAlgorAndCompare(fullfile(sInBaseDir, sMatFile), fullfile(sInBaseDir, sPosFile),...
%                                         cDatasetInfo{1}, cDatasetInfo{2}, vNum, cDatasetInfo{3}, cDatasetInfo{4}, cDatasetInfo{5}, cDatasetInfo{6}, cDatasetInfo{7},...
%                                             '', '', ccAlgors(c), 'sigmoidWeight', vSigmoidWeight(w), 'sigmoidMean', vSigmoidGamma(g), varargin{:});
%                 
%                                 assert(size(vBestObjVal,1) == 1);
% 
% 
% 
%                                 % store results
%                                 mResults( (r-1)*runsPerDataset + i, :) = [vBestObjVal(1), vComparisonDistNmi(1), vComparisonDistAmi(1), vComparisonDistAri(1), vRunningTime(1), vTotalIterTime(1), vDegenSolution(1)];
%                             end % end of loop through runsPerDataset
% 
%                         end % end of loop through numRum
%                 
%                         % compute average
%  
%                         vRows = find(mResults(:,degenCol) == 0);
%                         vAvg = mean(mResults(vRows,vRestCol),1);
%                         if length(vRows) <= 1
%                             vStd = zeros(1, length(vRestCol));
%                         else
%                             vStd = std(mResults(vRows,vRestCol),1);
%                         end
% %                         fprintf( '%d, %d, %.5f, .2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %d\n', vNum, p, vSigmoidWeight(w), vSigmoidGamma(g), vStats(1:5), vAvg, vStd, sum(mResults(:,degenCol) == 0));
%                         fprintf(fOutResults, '%d, %d, %.5f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %d\n', vNum, p, vSigmoidWeight(w), vSigmoidGamma(g), vStats(1:5), vAvg, vStd, sum(mResults(:,degenCol) == 0));                          
% 
%                     end % end of loop through sigmoidGamma
%                 end % end of loop through sigmoidWeight
%             
%             end % loop through number of positions
%                 
%         end % loop through density stats
%         
%         fclose(fOutResults);
%         % make sure all are closed
%         fclose('all');
%     
%     end % loop thru algorithms
    
    
end % end of function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bestObjVal, bestImageDis, comparisonDistNmi, comparisonDistAmi, comparisonDistAri, comparisonDistFuzzyAri, comparisonDistFuzzyJac, runningTime, totalIterNum, degenerate] =...
    runAndCompare(mAdj, mVertPosMapAct, mImageAct, cAlgor, varargin)
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

    
    % comparison measures between ground truth and obtained positions
    cfValidMeasFunc = {...
%         ClusterComparisonWrapper(@(a,b)(bmOverlapPosDist(a, mAdj, b, mAdj, 'nmi', 'discreteMemb', true)), mVertPosMapAct),...
%         ClusterComparisonWrapper(@(a,b)(bmOverlapPosDist(a, mAdj, b, mAdj, 'ami', 'discreteMemb', true)), mVertPosMapAct),...
%         ClusterComparisonWrapper(@(a,b)(bmOverlapPosDist(a, mAdj, b, mAdj, 'rand', 'discreteMemb', true)), mVertPosMapAct),...
    };    
    
    

    fImageDistanceFunc = EucImageDistance(mImageAct, mVertPosMapAct);
    fDistanceFunc = EucTriFact;
    fKKTDistanceFunc = EucKKTResidual(fDistanceFunc);

    display(cAlgor{2});
        
    if length(cAlgor) > 11
        cNewVarargin = cat(2, varargin, cAlgor{12:end});
    else
        cNewVarargin = varargin;
    end

        
    [mBestImage, mBestMembership, bestObjVal, ~, ~, ~, ~, ~, ~, runningTime, totalIterNum] =...
        runBMAlgor(mAdj, cAlgor{1:11}, fImageDistanceFunc, cfValidMeasFunc, cNewVarargin{:});
        
    %     assert(size(mBestMembership,2) == posNum);

        % check if the number of solutions obtained is less than expected.  If so,
        % we ignore that result and log that the number of clusters was less than
        % actipicated
    degenerate = 0;
    vColSum = sum(mBestMembership, 1);
    if sum(vColSum == 0) > 0
        warning('Algorithm %s has produced a degenerate solution with less than the number of specified clusters', cAlgor{2});
        degenerate = 1;
    end


    % compare the positions     
    comparisonDistNmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'nmi', 'discreteMemb', true);
    comparisonDistAmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'ami', 'discreteMemb', true);
    comparisonDistAri = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'rand', 'discreteMemb', true);
    comparisonDistFuzzyAri = fuzzyComparison(mVertPosMapAct, mBestMembership, 'ari');
    comparisonDistFuzzyJac = fuzzyComparison(mVertPosMapAct, mBestMembership, 'jaccard');

    bestImageDis = fImageDistanceFunc.distance(mBestImage, mBestMembership);
    kktResidual = fKKTDistanceFunc.residual(mAdj, mBestImage, mBestMembership);

end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% TODO: UPdate %
function loopThruPositions(cfOutResults, startPosNum, endPosNum, ccAlgors, sFileWildcard,...
    sInBaseDir, sMatSuffix, sPartSuffix, numRun, runsPerDataset,...
        bPosHeaderAct, bListFormatAct, vNum, bMatlabFormat, bSparse, bAddOneToIndex, bAddOneToPosVerts, bGraphHeader)

        % loop through the number of positions
        for p = startPosNum : endPosNum
            
            % number of results per run per algor
            resultNum = 10;
            cmResults = cell(size(ccAlgors,2),1);
            for c = 1 : size(ccAlgors,2)
                cmResults{c} = zeros(numRun, resultNum);
                ccAlgors{c}{3} = p;
            end
            
            % loop through the runs
            for r = 1 : numRun          
                
                % get the list of snapshot/matrices filenames
                sMatFileWildcard = sprintf('%s*%s', sFileWildcard, sMatSuffix);
                tsMatFilenames = dir(fullfile(sInBaseDir, sMatFileWildcard));
                
                % should only be one file
                assert(size(tsMatFilenames,1) == 1);

                sMatFile = tsMatFilenames(1).name;
                % get the associated position and image files
                sPosFile = strrep(sMatFile, sMatSuffix, sPartSuffix);
%                 sImageFile = strrep(sMatFile, sMatSuffix, sImageSuffix);
                
                for i = 1 : runsPerDataset

                    % no graph or comparison out file
                    [vBestObjVal, vBestMatApproxVal, vIdealEucObjVal, vIdealManObjVal, vComparisonDistVi, vComparisonDistAmi, vComparisonDistNmi, vRunningTime, vTotalIterTime, vDegenSolution] = ...
                            runBMAlgorAndCompare(fullfile(sInBaseDir, sMatFile), fullfile(sInBaseDir, sPosFile),...
                                bPosHeaderAct, bListFormatAct, vNum, bMatlabFormat, bSparse, bAddOneToIndex, bAddOneToPosVerts, bGraphHeader,...
                                    '', '', ccAlgors, varargin);
                
                    assert(size(vBestObjVal,1) == size(ccAlgors,2));
                
                    % store results
                    for c = 1 : size(ccAlgors,2)
                        cmResults{c}( (r-1)*runsPerDataset + i,:) = [vBestObjVal(c), vBestMatApproxVal(c), vIdealEucObjVal(c), vIdealManObjVal(c), vComparisonDistVi(c), vComparisonDistAmi(c), vComparisonDistNmi(c), vRunningTime(c), vTotalIterTime(c), vDegenSolution(c)];
                    end
                end

            end % end of loop through numRum
            
            % compute average
            for c = 1 : size(ccAlgors,2)
                % get all relevant columns that have non-degenerate solutions
                degenCol = 10;
                vRestCol = [1:9];
                vRows = find(cmResults{c}(:,degenCol) == 0);
                vAvg = mean(cmResults{c}(vRows,vRestCol),1);
                vStd = std(cmResults{c}(vRows,vRestCol),1);
                fprintf( '%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %d\n', vAvg, vStd, sum(cmResults{c}(:,degenCol) == 0));
                fprintf(cfOutResults{c}, '%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %d\n', vAvg, vStd, sum(cmResults{c}(:,degenCol) == 0));  
            end
            
            
        end % loop through number of positions

end

