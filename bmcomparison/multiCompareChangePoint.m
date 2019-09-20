function [mDistances, cChangePoints] = multiCompareChangePoint(sBaseDir, sMatWildcard, sSubseqFile, bSparse,...
    m, n, sOutFile, thres, bAddOne, sChangePointAlgor, bExit)
%
% Compares the conseqative blockmodels in the sequence of cPartitions to
% compare the performance of the various algorithms for different change
% point algorithms and parameters
%

    % load the sequence and partitions from file
    [cSnapshots, cPartitions, aPartInfo] = loadMatpart(sBaseDir, sMatWildcard, sSubseqFile, bSparse, m, n, bAddOne);


    % compute the different distances between adjacent blockmodels
    [mDistances] = multiCompareEnron(cSnapshots, cPartitions, aPartInfo);


    % compute change point
    cChangePoints = cell(size(mDistances,1),1);
    switch sChangePointAlgor
        case 'control'
            mu = mean(mDistances,2);
            % normalise by n
            stddev = std(mDistances, 0, 2);
            
            for r = 1 : size(mDistances,1)
                I = find(mDistances(r,:) > mu(r) + thres * stddev(r));
                cChangePoints{r} = I;
            end
            
        case 'cusum'
        
        otherwise
            warning('Unknown sChangePointAlgor = %s specified', sChangePointAlgor);
            return;
    
    end % end of switch
    
    % write out list of changepoints

    % open approprite file for writing
    fOut = fopen(sOutFile,'w');
    
    for l = 1 : size(cChangePoints)
        vChangePoints = cChangePoints{l};
        if size(vChangePoints,2) > 0
            for e = 1 : (size(vChangePoints,2)-1)
                fprintf(fOut, '%d,', vChangePoints(e));  
            end
            % print out last value
            fprintf(fOut, '%d\n', vChangePoints(size(vChangePoints,2)));
        end
    end
        
    fclose(fOut);
    
    if bExit
        exit;
    end
    
end % end of function