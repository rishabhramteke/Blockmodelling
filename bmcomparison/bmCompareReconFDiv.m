function [dist] = bmCompareReconFDiv(mMembership1, mMembership2, mImageMat1, mImageMat2, sDiv)
%
% mMembership1 - matrix of memberships of each vertex
% mMembership2 - matrix of memberships of each vertex
% mImageMat1 - image matrix
% mImageMat2 - image matrix
%
% Compares the blockmodels represented by cPosition1 and cPosition2,
% f-div reconstruction error.
%
% @author: Jeffrey Chan, 2013
%

    mAdj1 = mMembership1 * mImageMat1 * mMembership1';
    mAdj2 = mMembership2 * mImageMat2 * mMembership2';
    
    switch sDiv
        case 'kl'
            mAdj1Alt = mAdj1 + eps;
            mAdj2Alt = mAdj2 + eps;
            mF = klDiv(mAdj1Alt, mAdj2Alt);
            dist = sum(sum(mAdj1Alt .* mF));
        case 'hellinger'
            dist = 2 * sum(sum((sqrt(mAdj1) - sqrt(mAdj2)).^2));
        case 'pearson'
            dist = 0.5 * sum(sum(((mAdj1 - mAdj2).^2) ./ (mAdj2 + eps)));
        otherwise
            warning('Unknown');
    end
    
end % end of function


function [mF] = klDiv(mX, mY)
%
% KL distance between matrices.
%

    mZ = mY ./ mX;
    
    mF = mZ - log(mZ) - 1;

end % end of function