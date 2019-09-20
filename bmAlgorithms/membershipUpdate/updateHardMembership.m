function [mMembership] = updateHardMembership(mAdj, mImage, mMembership)
%
% Updates the hard membership, using a greedy, vertex by vertex assignment.
% Ordering of update does not affect final result.
%
% Batch update.
%
% Modified from Lo Bong's algorithm.
%

% To optimise (we can easily derive an update rule when a membership is
% moved to another)

% precompute some matrices
mMI = mMembership * mImage;
mIM = mImage * mMembership';
mExistApprox = mAdj - mMI * mMembership';

mOrigObjVal = trace(mExistApprox * mExistApprox');

% initialise two sparce vectors to quickly update memberships
mDelta = sparse(zeros(size(mMembership,1), size(mMembership,2)));

% initialise the new membership matrix
%mNewMembership = sparse(zeros(size(mMembership,1), size(mMembership,2)));


% loop through each vertex/row and assign to best position
for v = 1 : size(mMembership, 1)
%     disp(sprintf('vertex %d', v));
        
    % store objective values
    vObjVals = zeros(1,size(mMembership,2));

    
    % reset positions for mMembership
    %mMembership(v,:) = zeros(1,size(mMembership,2));
    % find which position v belongs to currently
    currP = find(mMembership(v,:) > 0);
    mDelta(v, currP) = -1;
    for p = 1 : size(mMembership,2)
        if (p == currP)
            vObjVals(p) = mOrigObjVal;
        else
            %mMembership(v,p) = 1;
            %vObjVals(p) = objective(mAdj, mImage, mMembership);
            %mMembership(v,p) = 0;
            mDelta(v, p) = 1;
            mDeltaTrans = mDelta';
            mNewApprox = mExistApprox + mMI * mDeltaTrans + mDelta * mIM - mDelta * mImage * mDeltaTrans;
            vObjVals(p) = trace(mNewApprox * mNewApprox');
            mDelta(v, p) = 0;
        end
    end % end of inner for
    mDelta(v, currP) = 0;
    
    % find minimum value
    minObjVal = min(vObjVals);
    vMinI = find(vObjVals == minObjVal);

    % only change assignment if minI != currP.
    if isempty(find(vMinI == currP, 1))
        % pick one of the minimum
        minP = vMinI(randsample(size(vMinI, 2), 1));
        % update mMembership
        mMembership(v, currP) = 0;
        mMembership(v,minP) = 1;
    end
   
end % end of outer for

end % end of function