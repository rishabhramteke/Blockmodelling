function [cmAlgorComparisonHypot] = evalFactoriseAlgors(cData, numRun, sOutBaseDir, sOutPrefix, varargin)
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

    % algorithm specifications (set default of one position, until they are set
    
    cccAlgors = {...
    
%             {'irm', 'irm', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'irm', bDiscretiseMembership, bColNorm},...
% %             {'bkn', 'bkn', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'bkn', bDiscretiseMembership, bColNorm},...        
%         {'s', {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCIncr', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'coordDescBMM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCIncr', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...        
%         {'s', {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'hardCIncrAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'hardCIncrAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'hardCIncrAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'coordDescBMM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'hardCIncrAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCPenaltyRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCPenaltyRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCPenaltyRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'hardCPenaltyRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'hardCPenaltyRowConstraintAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'hardCPenaltyRowConstraintAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'hardCPenaltyRowConstraintAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
%         {'s', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'bmEuclideanAdjEqual', 'hardCPenaltyRowConstraintAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...
          {'h', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runiteration', 100},   {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCPenaltyRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runiteration', 100, 'initPenaltyWeight', 1, 'penaltyMult', 10}},...                        
%         {'s', {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', 'randomInit', 'randomInit', true, bColNorm, 'runIterationNum', 100}},...        
%         {'s', {'matApprox', 'bolongM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', 'randomHardInit', 'randomHardInit', bDiscretiseMembership, bColNorm, 'runiteration', 400}},...
%         {'s', {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomHardInit', 'randomHardInit', true, bColNorm, 'runiteration', 400}},...
%         {'s', {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'hardCIncr', 'randomHardInit', 'randomHardInit', bDiscretiseMembership, bColNorm, 'runiteration', 400}},...
%         {'s', {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultOrthogonal', 'randomHardInit', 'randomHardInit', true, bColNorm, 'runiteration', 400}},...
%         {'s', {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', true, bColNorm, 'runiteration', 400}},...
%         {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'rowNormaliseAdhoc', true},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'rowNormaliseAdhoc', true},...          
%         {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'l1MemReg', 10},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'l1MemReg', 10},...
%         {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'l1ImageReg', 10},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'l1ImageReg', 10},...
%         {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'l1MemReg', 10, 'liImageReg', 10},...
% %         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'l1MemReg', 10, 'liImageReg', 10},...     
%         {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...
%         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...               
% %         {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', true, bColNorm},...
%             {'s', {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm}},...        
%             {'s', {'matApprox', 'multReparaM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0}},...
%             {'s', {'matApprox', 'multReparaM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0, 'l1MemReg', 100, 'l1ImageReg', 10}},...
%             {'s', {'matApprox', 'multLangM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0, 'l1MemReg', 100}},...            
%             {'s', {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.0}},... 
%             {'s', {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.5}},...                    
%                 {'s', {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 400}},...
%                 {'s', {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon+1, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 400, 'rowNormaliseAdhoc', true}},...
%                 {'s', {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 400}},...
%                 {'s', {'matApprox', 'dingM', initPosNum, convEpsilon+1, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 400, 'rowNormaliseAdhoc', true}},...                
%                 {'s', {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.0,  'runIterationNum', 100}},...             
%                 {'s', {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraintRepara', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.0,  'runIterationNum', 100}},... 
%                 {'s', {'matApprox', 'multLangM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0}},...
%                 {'s', {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'regWeight', 0, 'runIterationNum', 100, 'innerCollectPerIterStats', false}},...                        
%                 {'s', {'matApprox', 'coordDescExactBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclideanAdjEqual', 'softCMultHardRowConstraintAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.0,  'runIterationNum', 100}},...             
%                 {'s', {'matApprox', 'multLangM', initPosNum, convEpsilon+1, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0, 'l1MemReg', 100}},...                                            
%                 {'s', {'matApprox', 'multLangM', initPosNum, convEpsilon+1, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0, 'l1MemReg', 100}},...                            
%             {'s', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 10}},... 
%             {'s', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 10}},...
%             {'s', {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultHardRowConstraintRepara', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'runIterationNum', 100}},...              
%             {'h', {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'regWeight', 0, 'runIterationNum', 10},  {'matApprox', 'multLangM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0, 'runIterationNum', 1}},...                        
%             {'matApprox', 'multBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 1},...                    
%             {'matApprox', 'multReparaM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0},...
% {'bnmtf', 'bnmtfM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCBnmtf', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'regWeight', 1},...            
%             {'matApprox', 'coordDescBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.5},...                                
%         {'matApprox', 'bolongSoftMemM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultLong', 'randomInit', 'randomInit', true, bColNorm},...
% %         {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...        
%             {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'softCMultHardRowConstraint', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0},...
%             {'matApprox', 'multReparaM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0},...
%             {'matApprox', 'dingM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...        
%             {'matApprox', 'coordDescBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.5},...                    
%             {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.5},...                    
%             {'matApprox', 'coordDescExactM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 0.0},... 
% {'matApprox', 'multBMM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCMultDing', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm, 'sigmoidWeight', 1},...        
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'softCGradDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmEuclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...        
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCGradDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'bmCubicEuclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...                
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCGradDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclidean', 'softCCoordDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...
%         {'matApprox', 'projGradDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'softCCoordDescAdjEqual', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...
%         {'matApprox', 'coordDescM', initPosNum, convEpsilon, algorRunNum, 'euclideanAdjEqual', 'softCGradDesc', 'randomInit', 'randomInit', bDiscretiseMembership, bColNorm},...
        };
    

    % initialisations (same for all algorithms)
    sImageInit = 'randomInit';
    sMemInit = 'randomRowSumInit';    

    [cmAlgorComparisonHypot] = performRuns(cccAlgors, cData, numRun, sImageInit, sMemInit, sOutBaseDir, sOutPrefix, varargin{:});
    
    

    
    % close the files
%     for c = 1 : size(cData,2)
%         fclose(cfOutResults{c});
%     end
    
        
end % end of function



function [cmAlgorComparisonHypot] = performRuns(cccAlgors, cData, numRun, sImageInit, sMemInit, sOutBaseDir, sOutPrefix, varargin)


    % parse arguments
    inParser = inputParser;
    inParser.KeepUnmatched = true;    

    % default is to collect per iteration stats
    addParameter(inParser, 'collectPerIterStats', true);
    addParameter(inParser, 'hypothesisTesting', false);
    addParameter(inParser, 'maxStatsAcrossRuns', false);
    
    parse(inParser, varargin{:});
    
    bCollectPerIterStats = inParser.Results.collectPerIterStats;   
    bHypotTesting = inParser.Results.hypothesisTesting; 
    bMax = inParser.Results.maxStatsAcrossRuns;    
    

       
    
    
%     % initialise data structures
%     if bCollectPerIterStats
%         vBestObjVal = zeros(size(cccAlgors,2), 1);
%         vComparisonDistAmi = zeros(size(cccAlgors,2), 1);
%         vComparisonDistNmi = zeros(size(cccAlgors,2), 1);
%         vComparisonDistAri = zeros(size(cccAlgors,2), 1);
%         vRunningTime = zeros(size(cccAlgors,2), 1);
%         cmImage = cell(size(cccAlgors,2), 1);
%         cmMembership = cell(size(cccAlgors,2), 1);
%     end
    
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
                case 'randomRowSumInit'
                    fMemInitFunc = InitRowSumMem();
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
    
    
%     % construct outfiles
%     if ~bCollectPerIterStats
%         % one file per dataset
%         cfResultsFiles = cell(1, size(cccAlgors,2));
%         for c = 1 : size(cccAlgors,2)
%             % TODO: no hybrid output!
%             cAlgors = cccAlgors{c}{2};
%             sOutResultFile = sprintf('%s_%s_%s_%s_%.2f_%d%s', sOutPrefix, cAlgors{2}, cAlgors{7}, cAlgors{6}, cAlgors{4}, cAlgors{5},...
%                 '.results.csv');
%             cfResultsFiles{c} = fopen(fullfile(sOutBaseDir, sOutResultFile), 'a'); 
%         end    
%     end


    % PARAMETERS
    % number of results per run per algor
    resultNum = 11;
    % columns to get do hypthosise testing on
    vHypotCol = [1,4,5,9];


    cmAlgorComparisonHypot = {};
    % initialise datastructure
    if bHypotTesting && ~bCollectPerIterStats
        cmAlgorComparisonHypot = cell(1, length(vHypotCol));
        % algor by algor by dataset
        for col = 1 : length(vHypotCol)
            cmAlgorComparisonHypot{col} = zeros(size(cccAlgors,2),size(cccAlgors,2), size(cData,2));
        end  
    end

        
    
    
    % loop thru the datasets
    for d = 1 : size(cData,2)
        mAdj = cData{d}{1};
        mVertPosMapAct = cData{d}{2};
        mImageAct = cData{d}{3};
        sDataName = cData{d}{4};
        numPos = size(mVertPosMapAct,2);
        
        display(sDataName);
        
        % construct file 
        if ~bCollectPerIterStats
            sOutResultFile = sprintf('%s_%s%s', sOutPrefix, sDataName, '.results.csv');    
            fResultsFile =  fopen(fullfile(sOutBaseDir, sOutResultFile), 'a');
        end
        
        
        
        cmResults = cell(size(cccAlgors,2),1);
        for c = 1 : size(cccAlgors,2)
            cmResults{c} = zeros(numRun, resultNum);
            % update number of positions in the cccAlgors
            for b = 2 : size(cccAlgors{c},2)
                cccAlgors{c}{b}{3} = numPos;
            end
        end
        
        % comparison measures between ground truth and obtained positions
        if bCollectPerIterStats
            cfValidMeasFunc = {...
                ClusterComparisonWrapper(@(a,b)(fuzzyComparison(a,b,'rand')), mVertPosMapAct),...
                ClusterComparisonWrapper(@(a,b)(fuzzyComparison(a,b,'ari')), mVertPosMapAct),...
                ClusterComparisonWrapper(@(a,b)(bmOverlapPosDist(a, mAdj, b, mAdj, 'nmi', 'discreteMemb', true)), mVertPosMapAct),...
                ClusterComparisonWrapper(@(a,b)(bmOverlapPosDist(a, mAdj, b, mAdj, 'ami', 'discreteMemb', true)), mVertPosMapAct),...
            };
        else
            cfValidMeasFunc = {};
        end
        
 
%         fImageDistanceFunc = EucImageDistance(zeros(numPos, numPos), mVertPosMapAct); 
        
        fImageDistanceFunc = ImagePurtity();        
        fDistanceFunc = EucTriFact;
        fKKTDistanceFunc = EucKKTResidual(fDistanceFunc);    
        % modularity measures
%         fHardModularity = Modularity(mAdj, true);
%         fSoftModularity = Modularity(mAdj, false);
            
        
        % loop through the runs
        for r = 1 : numRun          
            % compare
%             [cmImage, cmMembership, vBestObjVal, vComparisonDistAmi, vComparisonDistNmi, vComparisonDistAri, vRunningTime] = ...
%                  runAndCompare(mAdj, mVertPosMapAct, mImageAct, sDataName, sImageInit, sMemInit, cfValidMeasFunc, cccAlgors, varargin{:});

            % initialise (this only tested for factorisation based algorithms)
            mInitMembership = fMemInitFunc.initMembership(mAdj, size(mVertPosMapAct,2));

            % initialise mImage
            mInitImage = fImageInitFunc.initImage(mAdj, mInitMembership);           

            algorRunNum = 1;
            % run each algorithm
            for a = 1 : size(cccAlgors,2)
                % determine if single or hybrid algorithm
                sAlgorType = cccAlgors{a}{1};
                switch sAlgorType
                    % single algorithm option
                    case 's'
                        ccAlgors = cccAlgors{a}{2};
%                         display(ccAlgors{2});

                        if length(ccAlgors) > 11
                            cNewVarargin = cat(2, varargin, ccAlgors{12:end});
                        else
                            cNewVarargin = varargin;
                        end

            
                        [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, runningTime, totalIterNum] =...
                             runBMAlgor(mAdj, ccAlgors{1}, ccAlgors{2}, ccAlgors{3}, ccAlgors{4}, algorRunNum, ccAlgors{6}, ccAlgors{7}, ccAlgors{8}, ccAlgors{9}, ccAlgors{10}, ccAlgors{11},...
                                 fImageDistanceFunc, cfValidMeasFunc, 'initialMem', mInitMembership, 'initialImage', mInitImage, 'collectPerIterStats', true, cNewVarargin{:});        


                    % hybrid algorithm option
                    case 'h'
                        ccAlgors1 = cccAlgors{a}{2};
                        ccAlgors2 = cccAlgors{a}{3};

                        if length(ccAlgors1) > 11
                            cNewVarargin1 = cat(2, varargin, {'initialMem', mInitMembership, 'initialImage', mInitImage, 'collectPerIterStats', true}, ccAlgors1{12:end});
                        else
                            cNewVarargin1 = varargin;
                        end     

                        if length(ccAlgors2) > 11
                            cNewVarargin2 = cat(2, varargin, {'collectPerIterStats', true}, ccAlgors2{12:end});
                        else
                            cNewVarargin2 = varargin;
                        end                           

                        [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, runningTime, totalIterNum] =...
                             binaryMembershipHyrbid(mAdj, ccAlgors1{1}, ccAlgors2{1}, ccAlgors1{2}, ccAlgors2{2}, ccAlgors1{3}, ccAlgors1{4}, ccAlgors1{4},...
                                 algorRunNum, ccAlgors1{6}, ccAlgors2{6}, ccAlgors1{7}, ccAlgors2{7}, ccAlgors1{8}, ccAlgors1{9}, ccAlgors1{10}, ccAlgors2{10},...
                                     ccAlgors1{11}, ccAlgors2{11}, fImageDistanceFunc, cfValidMeasFunc, cNewVarargin1, cNewVarargin2);                            
                    otherwise
                        error('evalFactoriseAlgors:runAndCompare', '%s is a unknown algorithm type in %dth specification for ccAlgors.\n',...
                            sAlgorType, a);
                end

                if sum(sum(~isreal(mBestMembership)))
%                     comparisonDistFuzzyAri = -1;
                    comparisonDistFuzzyJac = 0;
                    kktResidual = 10^6;
                else
%                     comparisonDistFuzzyAri = fuzzyComparison(mVertPosMapAct, mBestMembership, 'ari');
                    comparisonDistFuzzyJac = fuzzyComparison(mVertPosMapAct, mBestMembership, 'jaccard');
                    kktResidual = fKKTDistanceFunc.residual(mAdj, mBestImage, mBestMembership);
                end
                
                comparisonDistAmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'ami', 'discreteMemb', true);
                comparisonDistNmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'nmi', 'discreteMemb', true);
%                 comparisonDistAri = 1-bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'rand', 'discreteMemb', true);
                

%                 bestImageDis = fImageDistanceFunc.distance(mBestImage, mBestMembership);
%                 bestImageDis = 0;


                % modularity (hard and directed)
%                 bestHardModularity = fHardModularity.modularity(mBestMembership);
                % soft modularity (directed)
%                 bestSoftModularity = fSoftModularity.modularity(mBestMembership);
                
                % seqiBloc encoding cost
                [bestEncodingCost, bestAdjEncodingCost] = encodingCost(mAdj, mBestMembership);
                
                % number of zro rows
                zeroRowNum = length(find(sum(mBestMembership,2) == 0));
                

                % store results
%                 vImageHist = histc(full(cmImage{c}(:)), 0:0.1:0.1*(imageBinNum-1));
%                 vPosHist = histc(full(cmMembership{c}(:)), 0:0.1:0.1*(posBinNum-1));
                
                cmResults{a}(r,:) = [bestObjVal,  kktResidual, comparisonDistAmi, comparisonDistNmi, comparisonDistFuzzyJac,...
                    sum(sum(mBestMembership)), sum(sum(mBestImage)), bestEncodingCost, bestAdjEncodingCost, zeroRowNum, runningTime];

                
                
%                 if bCollectPerIterStats
%                     vBestObjVal(a) = bestObjVal;
%                     vComparisonDistAmi(a) = comparisonDistAmi;
%                     vComparisonDistNmi(a) = comparisonDistNmi;
%                     vComparisonDistAri(a) = comparisonDistAri;
%                     vRunningTime(a) = runningTime;
%                     cmImage{a} = mBestImage;
%                     cmMembership{a} = mBestMembership;
%                 end

                % display the various stats about each run of each algorithm 
                % (should have same number of positions whether single or hybrid, so
                % cccAlgors{a}{2}{3} = posNum).
                if bCollectPerIterStats
                    plotInfo(sDataName, cccAlgors{a}{2}{3}, mAdj, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmMembership, ccmImage);
                    mBestMembership
                    mBestImage
                    runningTime
                    sum(mBestMembership, 2)
                end

            end  % end of for loop through algorithms              
                
%             if bCollectPerIterStats
%                 assert(size(vBestObjVal,1) == size(cccAlgors,2));
%             end
    

        end % end of loop through numRum
        
%         % sort by nmi
%         for c = 1 : size(cccAlgors,2)
%             [~, vSortedIndex] = sort(cmResults{c}(:,3), 'descend');
%             cmResults{c} = cmResults{c}(vSortedIndex, :);
%         end



        % hypthoesis testing
        if bHypotTesting && bCollectPerIterStats
             % -1 is loss to our compared algorithm, 0 is draw, 1 is win
             
             
        end

        % average and output
        if ~bCollectPerIterStats
            for c = 1 : size(cccAlgors,2)
                nmiCol = 4;
                if bMax
                    % find max
                    [~, maxIndex] = max(cmResults{c}(:,nmiCol));
                    vAvg = cmResults{c}(maxIndex,:);
                else
                    % take average across runs
                    vAvg = mean(cmResults{c},1);
                end
                
                if numRun <= 1
                    vStd = zeros(1, size(cmResults{c},2));
                else
                    vStd = std(cmResults{c},1);
                end
                cAlgors = cccAlgors{c}{2};
                fprintf('%s, %s, %s, %d, %.2f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f,  %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n', cAlgors{2}, cAlgors{7}, cAlgors{6}, cAlgors{4}, cAlgors{5}, vAvg, vStd);                          
                fprintf(fResultsFile, '%s, %s, %s, %d, %.2f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n', cAlgors{2}, cAlgors{7}, cAlgors{6}, cAlgors{4}, cAlgors{5}, vAvg, vStd);                          


            end        
            
            drawAlpha = 0.1;
            oneSidedAlpha = 0.05;
            
            % hypthoesis testing
            if bHypotTesting
                for c1 = 1 : size(cccAlgors,2)
                    for c2= c1 : size(cccAlgors,2)
                        for hycol = 1 : length(vHypotCol)
                            currCol = vHypotCol(hycol);
%                             fprintf('algor %d %d %d\n', c1, c2, currCol);
                            vRowsToUse = 1:1:length(cmResults{c1}(:,currCol));
                            vRowToRemove = [];
                            if sum(isnan(cmResults{c1}(:,currCol))) > 0 || sum(isinf(cmResults{c1}(:,currCol))) > 0
                                % remove the row
                                if sum(isnan(cmResults{c1}(:,currCol))) > 0
                                    vI = find(isnan(cmResults{c1}(:,currCol)));
                                else
                                    vI = find(isinf(cmResults{c1}(:,currCol)));
                                end
                                vRowToRemove = [vRowToRemove vI];
%                                 cmResults{c1}(:,currCol)
                            end
                            if sum(isnan(cmResults{c2}(:,currCol))) > 0 || sum(isinf(cmResults{c2}(:,currCol))) > 0
                                % remove the row
                                if sum(isnan(cmResults{c2}(:,currCol))) > 0 
                                    vI = find(isnan(cmResults{c2}(:,currCol)));
                                else
                                    vI = find(isinf(cmResults{c2}(:,currCol)));
                                end
                                vRowToRemove = [vRowToRemove vI];
                            end
                            vRowsToUse(vRowToRemove) = [];
                            vResults1 = real(cmResults{c1}(vRowsToUse,currCol));
                            vResults2 = real(cmResults{c2}(vRowsToUse,currCol));
                            % test if draw first
                            [hypothesis, pValue] = ttest(vResults1, vResults2, 'Alpha', drawAlpha, 'Tail', 'both');
                            if hypothesis == 0
                                [hypothesis, pValue] = ttest(vResults1, vResults2, 'Alpha', oneSidedAlpha, 'Tail', 'right');
                                if hypothesis == 0
                                    [hypothesis, pValue] = ttest(vResults1, vResults2, 'Alpha', oneSidedAlpha, 'Tail', 'left');
                                    if hypothesis == 0
                                        % inconclusive
                                        finalHypotValue = 2;
                                    else
                                        % less than
                                        finalHypotValue = -1;
                                    end
                                else
                                    % greater than (for the first one)
                                    finalHypotValue = 1;
                                end
                            else
                                % draw!
                                finalHypotValue = 0;
                            end
                                
                            % -1 if mean of first set of results is less than the 2nd set, 0 is draw, 1
                            % is greater, 2 is inconclusive


                            cmAlgorComparisonHypot{hycol}(c1,c2,d) = finalHypotValue;
                        end
                    end
                end
            end            
        end
        
        if ~bCollectPerIterStats
            fclose(fResultsFile);
        end
        
        

    end % end of loop thru datasets


end % end of function


% function [cmImage, cmMembership, vBestObjVal, vComparisonDistAmi, vComparisonDistNmi, vComparisonDistAri, vRunningTime] =...
%     runAndCompare(mAdj, mVertPosMapAct, mImageAct, sDataName, sImageInit, sMemInit, cfValidMeasFunc, cccAlgors, varargin)
% 
%     % check if we need to plot and store the iteration stats
%     
%     % parse arguments
%     inParser = inputParser;
%     inParser.KeepUnmatched = true;    
% 
%     % default is to collect per iteration stats
%     addParameter(inParser, 'collectPerIterStats', true);
%     
%     parse(inParser, varargin{:});
%     
%     bCollectPerIterStats = inParser.Results.collectPerIterStats;   
%     
%     
%     % initialise data structures
%     if bCollectPerIterStats
%         vBestObjVal = zeros(size(cccAlgors,2), 1);
%         vComparisonDistAmi = zeros(size(cccAlgors,2), 1);
%         vComparisonDistNmi = zeros(size(cccAlgors,2), 1);
%         vComparisonDistAri = zeros(size(cccAlgors,2), 1);
%         vRunningTime = zeros(size(cccAlgors,2), 1);
%         cmImage = cell(size(cccAlgors,2), 1);
%         cmMembership = cell(size(cccAlgors,2), 1);
%     end
%     
%     % image initialisation
%     switch sImageInit
%         case 'randomInit'
%             fImageInitFunc = InitRandomImage(false);
%         case 'randomHardInit'
%             fImageInitFunc = InitRandomImage(true);
%         case 'kmeans'
%             fImageInitFunc = InitKmeansImage();
%         case 'svd'
%             fImageInitFunc = InitSvdImage();
%         otherwise 
%             error('runAndCompare:sImageInit', 'Invalid image initialisation option');
%     end % end of switch    
%     
%     
%     % membership initialisation
%     switch sMemInit
%         case 'randomInit'
%             fMemInitFunc = InitRandomMem(false);
%         case 'randomHardInit'
%             fMemInitFunc = InitRandomMem(true);
%         case 'kmeans'
%             fMemInitFunc = InitKmeansMem(true);
%         case 'svd'
%             fMemInitFunc = InitSvdMem(true);
%         case 'randomAcol'
%             fMemInitFunc = InitRandomAcolMem();
%         otherwise 
%             error('runAndCompare:sMemInit', 'Invalid membership initialisation option');
%     end % end of switch          
%     
%     
%     % initialise (this only tested for factorisation based algorithms)
%     mInitMembership = fMemInitFunc.initMembership(mAdj, size(mVertPosMapAct,2));
% 
%     % initialise mImage
%     mInitImage = fImageInitFunc.initImage(mAdj, mInitMembership);    
%     
%     fImageDistanceFunc = EucImageDistance(mImageAct, mVertPosMapAct);
%         
%     algorRunNum = 1;
%     % run each algorithm
%     for a = 1 : size(cccAlgors,2)
%             % determine if single or hybrid algorithm
%             sAlgorType = cccAlgors{a}{1};
%             switch sAlgorType
%                 % single algorithm option
%                 case 's'
%                     ccAlgors = cccAlgors{a}{2};
%                     display(ccAlgors{2});
%                      
%                     
%                     if length(ccAlgors) > 11
%                         cNewVarargin = cat(2, varargin, ccAlgors{12:end});
%                     else
%                         cNewVarargin = varargin;
%                     end
%                         
% %                         [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, runningTime, totalIterNum] =...
% %                             runBMAlgor(mAdj, ccAlgors{a}{1}, ccAlgors{a}{2}, ccAlgors{a}{3}, ccAlgors{a}{4}, algorRunNum, ccAlgors{a}{6}, ccAlgors{a}{7}, ccAlgors{a}{8}, ccAlgors{a}{9}, ccAlgors{a}{10}, ccAlgors{a}{11},...
% %                                 fImageDistanceFunc, cfValidMeasFunc, , 'initialMem', mInitMembership, 'initialImage', mInitImage, 'collectPerIterStats', true, varargin{:});        
% %                    
%                     [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, runningTime, totalIterNum] =...
%                          runBMAlgor(mAdj, ccAlgors{1}, ccAlgors{2}, ccAlgors{3}, ccAlgors{4}, algorRunNum, ccAlgors{6}, ccAlgors{7}, ccAlgors{8}, ccAlgors{9}, ccAlgors{10}, ccAlgors{11},...
%                              fImageDistanceFunc, cfValidMeasFunc, 'initialMem', mInitMembership, 'initialImage', mInitImage, 'collectPerIterStats', true, cNewVarargin{:});        
% 
%                            
%                 % hybrid algorithm option
%                 case 'h'
%                     
%                     ccAlgors1 = cccAlgors{a}{2};
%                     ccAlgors2 = cccAlgors{a}{3};
% 
%                     if length(ccAlgors1) > 11
%                         cNewVarargin1 = cat(2, varargin, {'initialMem', mInitMembership, 'initialImage', mInitImage, 'collectPerIterStats', true}, ccAlgors1{12:end});
%                     else
%                         cNewVarargin1 = varargin;
%                     end     
%                     
%                     if length(ccAlgors2) > 11
%                         cNewVarargin2 = cat(2, varargin, {'collectPerIterStats', true}, ccAlgors2{12:end});
%                     else
%                         cNewVarargin2 = varargin;
%                     end                           
%                     
%                     [mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, runningTime, totalIterNum] =...
%                         binaryMembershipHyrbid(mAdj, ccAlgors1{1}, ccAlgors2{1}, ccAlgors1{2}, ccAlgors2{2}, ccAlgors1{3}, ccAlgors1{4}, ccAlgors1{4},...
%                             algorRunNum, ccAlgors1{6}, ccAlgors2{6}, ccAlgors1{7}, ccAlgors2{7}, ccAlgors1{8}, ccAlgors1{9}, ccAlgors1{10}, ccAlgors2{10},...
%                                 ccAlgors1{11}, ccAlgors2{11}, fImageDistanceFunc, cfValidMeasFunc, cNewVarargin1, cNewVarargin2);                            
%                     
%                 otherwise
%                     error('evalFactoriseAlgors:runAndCompare', '%s is a unknown algorithm type in %dth specification for ccAlgors.\n',...
%                         sAlgorType, a);
%             end
%             
%             comparisonDistAmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'ami');
%             comparisonDistNmi = bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'nmi');
%             comparisonDistAri = 1-bmOverlapPosDist(mVertPosMapAct, mAdj, mBestMembership, mAdj, 'rand');
% 
%         if bCollectPerIterStats
%             vBestObjVal(a) = bestObjVal;
%             vComparisonDistAmi(a) = comparisonDistAmi;
%             vComparisonDistNmi(a) = comparisonDistNmi;
%             vComparisonDistAri(a) = comparisonDistAri;
%             vRunningTime(a) = runningTime;
%             cmImage{a} = mBestImage;
%             cmMembership{a} = mBestMembership;
%         end
%         
%         % display the various stats about each run of each algorithm 
%         % (should have same number of positions whether single or hybrid, so
%         % cccAlgors{a}{2}{3} = posNum).
%         if bCollectPerIterStats
%             plotInfo(sDataName, cccAlgors{a}{2}{3}, mAdj, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmMembership, ccmImage);
%             mBestMembership
%             mBestImage
%         end
%             
%     end  % end of for loop   
%         
% 
%         
% end % end of function


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
