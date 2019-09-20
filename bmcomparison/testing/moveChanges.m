function [cmResults] = moveChanges(mAdj, mInitMembership, percChange, numRuns, ccMeasures)
%
% @author: Jeffrey Chan, 2013

    % generate changes, up to percChange
    vertNum = size(mAdj,1);
    changeNum = max(floor(vertNum * percChange), 1);

    cmResults = cell(1, length(ccMeasures));
    for m = 1 : length(ccMeasures)
        cmResults{m} = zeros(numRuns, changeNum);
    end
    
    posNum = size(mInitMembership,2);
    cvPos = cell(1, posNum);
    for p = 1 : posNum
        vR = find(mInitMembership(:,p));
        cvPos{p} = vR;
    end
    

    % loop through runs
    for r = 1 : numRuns
        mMembership = mInitMembership;
        
        
        vPosSize = sum(mInitMembership, 1);
        
        vMoves = randperm(vertNum);
        c = 1;
        while c <= changeNum
            % check if we can move
            currPos = find(mMembership(vMoves(c),:));
            if vPosSize(currPos) > 1
                newPos = randi(posNum-1);
                if newPos >= currPos
                    newPos = newPos + 1;
                end
                
                mMembership(vMoves(c), currPos) = 0;
                mMembership(vMoves(c), newPos) = 1;
                
                % loop through the measures
                for m = 1 : length(ccMeasures)
                    switch ccMeasures{m}{2}
                        case 1
                            dist = bmCompareWeightedRecon(mInitMembership, mMembership, mAdj, mAdj, ccMeasures{m}{1});
                        case 2
                            dist = bmOverlapPosDist(mInitMembership, mAdj, mMembership, mAdj, ccMeasures{m}{1});
                        case 3
                            dist = bmOverlapDenDist(mInitMembership, mAdj, mMembership, mAdj, ccMeasures{m}{1});
                        otherwise
                            warning(sprintf('Unknown measure option %d for ccMeasure', ccMeasures{m}{2}));
                    end
                    cmResults{m}(r, c) = dist;
                end
                
                vPosSize(newPos) = vPosSize(newPos) + 1;
                vPosSize(currPos) = vPosSize(currPos) - 1;
                c = c + 1;
            else
                % remove the move we can't make
                vMoves(c) = [];
            end
            
        end
        
    end % end of for
    

end % end of fucntion