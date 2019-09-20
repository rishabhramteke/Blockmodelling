function [mMembership] = updateHardMembershipInc(mAdj, mImage, mMembership, vargin)
%
% Updates the hard membership, using a greedy, vertex by vertex assignment.
% Ordering of update does not affect final result.
%
% Incremental update.
%

% To optimise (we can easily derive an update rule when a membership is
% moved to another)


bHasAdj = false;

if len(vargin) > 1
    mAdjWeight = vargin{1};
    bHasAdj = true;
end


% precompute some matrices
mMI = mMembership * mImage;
mIM = mImage * mMembership';
if ~bHasAdj
    mExistApprox = mAdj - mMI * mMembership';
else
    mExistApprox = (mAdj - mMI * mMembership') .* mAdjWeight;
end




currObjVal = trace(mExistApprox * mExistApprox');

% initialise two sparce vectors to quickly update memberships
mDelta = sparse(zeros(size(mMembership,1), size(mMembership,2)));

% initialise the new membership matrix
%mNewMembership = sparse(zeros(size(mMembership,1), size(mMembership,2)));


% randomise order we evaluate the vertices
vVertOrder = randperm(size(mMembership,1));

% loop through each vertex/row and assign to best position
for i = 1 : size(vVertOrder, 2)
%     disp(sprintf('vertex %d', vVertOrder(i)));
        
    % store objective values
    vObjVals = zeros(1,size(mMembership,2));

    % reset positions for mMembership
    % mMembership(v,:) = zeros(1,size(mMembership,2));
    % find which position v belongs to currently
    currP = find(mMembership(vVertOrder(i),:) > 0);
%     assert(~isempty(currP));
    % if singleton, we don't try moving
    if sum(mMembership(:, currP)) <= 1
        continue;
    end

    mDelta(vVertOrder(i), currP) = 1;
    
    for p = 1 : size(mMembership,2)
        if (p == currP)
            vObjVals(p) = currObjVal;
        else
            %mMembership(v,p) = 1;
            %vObjVals(p) = objective(mAdj, mImage, mMembership);
            %mMembership(v,p) = 0;
            mDelta(vVertOrder(i), p) = -1;
            mDeltaTrans = mDelta';
            mNewApprox = mExistApprox + (mMI * mDeltaTrans + mDelta * mIM - mDelta * mImage * mDeltaTrans) .* mAdjWeight;
            vObjVals(p) = trace(mNewApprox * mNewApprox');
            mDelta(vVertOrder(i), p) = 0;
        end
    end % end of inner for
    mDelta(vVertOrder(i), currP) = 0;
    
    % find minimum value
    minObjVal = min(vObjVals);
    vMinI = find(vObjVals == minObjVal);
%     assert(size(vMinI,2) > 0);

    % only change assignment if minI != currP 
    if isempty(find(vMinI == currP, 1))
        % pick one of the minimum
        minP = vMinI(randsample(size(vMinI,2), 1));
        % update mMembership
        mMembership(vVertOrder(i), currP) = 0;
        mMembership(vVertOrder(i), minP) = 1;
        % update all the precalculated matrices
        mDelta(vVertOrder(i), minP) = -1;
        mDeltaTrans = mDelta';
        mExistApprox = mExistApprox + (mMI * mDeltaTrans + mDelta * mIM - mDelta * mImage * mDeltaTrans) .* mAdjWeight;
        mMI = mMI - mDelta * mImage;
        mIM = mIM - mImage * mDeltaTrans;
        % update objective value
        currObjVal = minObjVal;
    end
   
end % end of outer for

end % end of function