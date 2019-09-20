function [mApprox] = genBlockmodelApproxGraph(mAdj, mMembership, mImage)
%
% Generates a blockmodel approximated graph.
%

    % find the image matrix according to mAdj and mMembership
    posNum = size(mMembership,2);
    cvMembership = cell(1, posNum);
    for p = 1 : posNum
       cvMembership{p} = find(mMembership(:,p) > 0); 
    end
    
    mApprox = zeros(size(mAdj,1), size(mAdj,2));
    
    % expected (density) image matrix
    mExpImage = zeros(posNum, posNum);
    for r = 1 : posNum
        for c = 1 : posNum
%             cvMembership{r}
%             cvMembership{c}
%             mAdj(cvMembership{r}, cvMembership{c})
            mExpImage(r,c) = length(nonzeros(mAdj(cvMembership{r}, cvMembership{c}))) /...
                (length(cvMembership{r}) * length(cvMembership{c}));
%             mRandom = rand(length(cvMembership{r}), length(cvMembership{c}));
%             mEdges = mRandom < mExpImage(r,c);
            mEdges = ones(length(cvMembership{r}), length(cvMembership{c})) * mExpImage(r,c);
            mApprox(cvMembership{r}, cvMembership{c}) = mEdges;
        end
    end
    
    mExpImage
    mApprox = sparse(mApprox);


end % end of function