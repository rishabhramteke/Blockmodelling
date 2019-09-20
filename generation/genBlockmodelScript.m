function genBlockmodelScript(sFileoutPrefix, numGraphPerTest, vGraphSize, vPosNum,...
    sPosDist, sBlockmodelType, vMeanDense, vStdDense, vBackgroundProp, bShuffle, bFileOutput, bEdgeList, varargin)
%
% Generates a series of graphs that have the same sets of positions for each
% graph size, graph run and position number.  The densities of the generated, with background noise portion, is
% varied for the same set of positions.  
%
% See code for details about arguments.
%
% @author Jeffrey Chan
% @Date 2013
%

    vPosNumStart = vPosNum(1);
    vPosNumEnd = vPosNum(2);
    if vPosNumStart > vPosNumEnd
       error('Range of the number of positions is incorrect'); 
    end
    

    
    for g = 1 : length(vGraphSize)
        graphSize = vGraphSize(g);
        fprintf('Generating graph size %d.\n', graphSize);
        for i = 1 : numGraphPerTest
            fprintf('Generating graph set %d.\n', i);
            for posNum = vPosNumStart : vPosNumEnd
                fprintf('Generating posNum = %d.\n', posNum);
                % unweighted graph
                genUnweighted(sFileoutPrefix, sPosDist, sBlockmodelType,...
                    vMeanDense, vStdDense, vBackgroundProp, bShuffle, bFileOutput, bEdgeList, i, graphSize, posNum, varargin{:});            
            end
        end
        
    end % end of for loop through graph sizes
             
end % end of function



function genUnweighted(sFileoutPrefix, sPosDist, sBlockmodelType,...
            vMeanDense, vStdDense, vBackgroundProp, bShuffle, bFileOutput, bEdgeList, graphIndex, graphSize, posNum, varargin)
    %
    % Generate unweighted graphs.
    %
    
    
    % file extension type
    sAdjFileExtension = '.vel.mat.csv';
    sPosFileExtension = '.part.csv';
    sImageFileExtension = '.image.csv';
    
    % defaults
    meanSparse = 0.0;
    stdSparse = 0.0;
    
    
    if length(varargin) > 0
        [~, ~, ~, cmImageGraph, cmShuffledAdjMat, vShuffledVertPos] =...
            genBlockmodelWithBackground(sPosDist, posNum, graphSize, bShuffle, vBackgroundProp, sBlockmodelType, vMeanDense, vStdDense, varargin{1});
    else
        [~, ~, ~, cmImageGraph, cmShuffledAdjMat, vShuffledVertPos] =...
            genBlockmodelWithBackground(sPosDist, posNum, graphSize, bShuffle, vBackgroundProp, sBlockmodelType, vMeanDense, vStdDense);
    end

    assert(length(vMeanDense) == size(cmShuffledAdjMat, 2));
            
    % loop through each of the dense entries and construct and write out
    % the appropriate files
    for k = 1 : length(vMeanDense)
        for b = 1 : length(vBackgroundProp)
             sFilenameSuffix = sprintf('%d-%d-%d-%s-%s-%.2f-%.2f-%.2f-%.2f-%.2f',...
                 graphIndex, graphSize, posNum, sPosDist, sBlockmodelType, vBackgroundProp(b), vMeanDense(k), vStdDense(k), meanSparse, stdSparse);
             if length(varargin) > 0
                 sFilenameSuffix = strcat(sFilenameSuffix, sprintf('-opt-%.2f', varargin{1}));
             end
                
             sAdjFilename = strcat(sFileoutPrefix, '-', sFilenameSuffix, sAdjFileExtension);
             sPosFilename = strcat(sFileoutPrefix, '-', sFilenameSuffix, sPosFileExtension);
             sImageFilename = strcat(sFileoutPrefix, '-', sFilenameSuffix, sImageFileExtension);
                
             % if we want to write out to file
             % otherewise code doesn't generate any output
             if bFileOutput
                 if bEdgeList
                    [vR, vC] = find(cmShuffledAdjMat{b,k});
                    mComb = [vR-1 vC-1];
                    dlmwrite(sAdjFilename, mComb, 'precision', '%d');
                 else
                    dlmwrite(sAdjFilename, cmShuffledAdjMat{b,k}, 'precision', '%d');
                 end
    
                 % loop through the partitions and write out to file (note it appends)
                 for p = 1 : posNum
                     vIndices = find(vShuffledVertPos == p);
                     dlmwrite(sPosFilename, vIndices - 1, '-append'); 
                 end
    
                 dlmwrite(sImageFilename, cmImageGraph{b,k}, 'precision', '%.2f');
             end                
        end
    end    
end % end of function

