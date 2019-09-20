function [posDis] = bmOverlapPosDist(mMembership1, mAdj1, mMembership2, mAdj2, sPosDist, varargin)
%
% Computes the positional distances between clustering/bm 1 and
% 2.
%
% sPosDist - Position distance to use.
%

% mMembership1
% mMembership2

    % parse arguments
    inParser = inputParser;

    % add parameters and set default valuse
    addParameter(inParser, 'discreteMemb', false);

    parse(inParser, varargin{:});
    
    bDiscretiseMembership = inParser.Results.discreteMemb;



% convert to vector of membership
% TODO: assuming hard partitioning
vMembership1 = zeros(size(mMembership1,1), 1);
vMembership2 = zeros(size(mMembership2,1), 1);
for v = 1 : size(mMembership1,1)
    assert(sum(mMembership1(v,:)) > 0);
    if bDiscretiseMembership
        [~,pos] = max(mMembership1(v,:));
        vMembership1(v) = pos;
    else
        % TODO: this doesn't make sense to just use first non-zero entry without
        % checking that there is only one non-zero value
        vPos = find(mMembership1(v,:) > 0);
        vMembership1(v) = vPos(1);
    end
end
for v = 1 : size(mMembership2,1)
%     assert(sum(mMembership2(v,:)) > 0);
    if sum(mMembership2(v,:)) == 0
        % assign the opposite position to vMembership1(v)
        vPossPos = 1:size(mMembership2,2);
        % delete possible position that  corresponds to the one in vMembership1
        vPossPos(vMembership1(v)) = [];
        vMembership2(v) = randsample(vPossPos, 1);
       % assign random if this occurs and print warning
%        vMembership2(v) = randi(size(mMembership2,2), 1);
       warning('Vertex %d does not have a position membership.  Random membership assigned to it (that is not same as vMembership1).\n', v);
    else
        if bDiscretiseMembership
            [~,pos] = max(mMembership2(v,:));
            vMembership2(v) = pos;
        else
            % TODO: this doesn't make sense to just use first non-zero entry without
            % checking that there is only one non-zero value
            vPos = find(mMembership2(v,:) > 0);
            vMembership2(v) = vPos(1);
        end
    end
end



switch sPosDist
    case 'vi'
        posDis = varOfInfoCopy(vMembership1, vMembership2);
        % normalise
        posDis = posDis / log2(size(vMembership1,1));
    case 'ami'
        % adjusted for chance mi (Vinh's measure)
        posDis = ami(vMembership1, vMembership2);
    case 'nmiMustHaveSameNoOfPos'
        % normalised mutual informaiton, where there must be the same number of
        % positions.  This version is faster than nmi2().
        posDis = nmi(vMembership1, vMembership2);
    case 'nmi'
        % normalised mutual informaiton, where there is no restrictions on the
        % number of positions.
        posDis = nmi2(vMembership1, vMembership2);
    case 'rand'
        [AR,~,~,~] = randIndex(vMembership1, vMembership2);
        % we use adjusted rand (AR), and convert it to a
        % distance/dissimilarity
        posDis = 1-AR;
    % the following two aren't really position distances, but included
    % because we need to use the vMemberships
    case 'jaccard'
        [~, ~, Jac, ~] = valid_external(vMembership1, vMembership2);
        posDis = 1 - Jac;
    case 'esVi'
        % we want non-edit distance version
        posDis = edgeSetVarOfInfo(vMembership1, mAdj1, vMembership2, mAdj2, 0);
        % normalise
        posDis = posDis / (2*log2(size(vMembership1,1)));
    case 'esRand'
        [AR,~,~,~] = edgeSetRandIndex(vMembership1, mAdj1, vMembership2, mAdj2);
        posDis = 1 - AR;
    otherwise
        warning('Unknown sPosDist = %s specified', sPosDist);
        return;
end



end % end of function




