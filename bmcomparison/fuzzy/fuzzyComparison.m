function [val] = fuzzyComparison(mMembership1, mMembership2, sBaseMeasure, varargin)
%
% Fuzzy cluster comparison measure.
%
% @author: Jeffrey Chan, 2014
%


    inParser = inputParser;
    inParser.KeepUnmatched = true;

    % add parameters and set default valuse
    addParameter(inParser, 'normalise', false);

    parse(inParser, varargin{:});
    
    bNorm = inParser.Results.normalise;

    % if normalise, we make sure each row of both membership sum to 1
    % if any row has all zeros, we randomly assign a number that is different
    % from the other
    if bNorm
        mNewMembership2 = mMembership2;
        for v = 1 : size(mMembership2,1)
            if sum(mMembership2(v,:)) == 0
                % assign the opposite position to vMembership1(v)
                vPossPos = 1:size(mMembership2,2);
                % delete possible position that  corresponds to the one in vMembership1
                [~, otherPos] = max(mMembership1(v,:));
                vPossPos(otherPos) = [];
                mNewMembership2(v, randsample(vPossPos, 1)) = 1;
                warning('Vertex %d does not have a position membership.  Random membership assigned to it (that is not same as vMembership1).\n', v);
            else
                rowTotal = sum(mMembership2(v,:));
                mNewMembership2(v,:) = mNewMembership2(v,:) / rowTotal; 
            end
        end
    end
    
    
    objNum = size(mMembership1, 1);

    mContingency = mMembership1' * mMembership2;

    
    scalingFactor =  objNum/ sum(sum(mContingency));
    mContingency = mContingency * scalingFactor;
    
    % agreement and disagreement sums
    vRowSum = sum(mContingency,2);
    vColSum = sum(mContingency,1);
    totalSum = sum(sum(mContingency));
    totalSqSum = sumsqr(mContingency);
    rowSqSum = sumsqr(vRowSum);
    colSqSum = sumsqr(vColSum);
    
    a = 0.5 * (totalSqSum - totalSum);
    d = 0.5 * (objNum^2 + totalSqSum - rowSqSum - colSqSum);
    b = 0.5 * (colSqSum - totalSqSum);
    c = 0.5 * (rowSqSum - totalSqSum);
    
    switch sBaseMeasure
        case 'rand'
            val = (a + d) / (a + b + c + d);
        case 'ari'
            val = (a - ((a+c) * (a+b) / (a + b + c + d))) / (0.5 * ((a+c) + (a+b)) - ((a+c)*(a+b) / (a+b+c+d)));
        case 'jaccard'
            val = a / (a + b + c);
        
        otherwise
            warning('fuzzyComparison:blah', 'sBaseMeasure %s is unknown.', sBaseMeasure);
    end % end of switch


end % end of function