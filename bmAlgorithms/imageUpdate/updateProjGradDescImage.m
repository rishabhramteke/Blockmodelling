function [mImageNew] = updateProjGradDescImage(mAdj, mImage, mMembership, fObjective)
%
% Updates mImage according to a projected gradient descent update rule.
%
%
% @author: Jeffrey Chan, 2013
%
    
    % start off roughly in the scale of [0,1] for M
    alpha = (size(mImage,1) / size(mAdj,1))^2;
    beta = 0.1;
    sigma = 0.01;
    
    mGrad = fObjective.gradient(mAdj, mImage, mMembership);
    mImageNew = posProjection(mImage + alpha * mGrad);
    currDis = fObjective.distance(mAdj, mImage, mMembership);
    
%     if fObjective.distance(mAdj, mImageNew, mMembership) - currDis > sigma * mGrad' * (mImageNew - mImage)
    if fObjective.distance(mAdj, mImageNew, mMembership) - currDis > sigma * sum(sum(mGrad .* (mImageNew - mImage)))
        % compute new alpha and new image and gradient
        alpha = alpha * beta;
        mImageNew = posProjection(mImage + alpha * mGrad);
        % while difference in magnitude isn't a "minimum", decrease alpha
        while fObjective.distance(mAdj, mImageNew, mMembership) - currDis > sigma * sum(sum(mGrad .* (mImageNew - mImage)))
%             disp(sprintf('alpha = %f', alpha));
            % compute new alpha and new image and gradient
            alpha = alpha * beta;
            mImageNew = posProjection(mImage + alpha * mGrad);
        end
    else
        % compute new alpha and new image and gradient
        alpha = alpha / beta;
        imageDiffLimit = 10^-9;
        mImageNew = posProjection(mImage + alpha * mGrad);
        while abs(sum(sum((mImageNew - mImage)))) > imageDiffLimit && (fObjective.distance(mAdj, mImageNew, mMembership) - currDis <= sigma * sum(sum(mGrad .* (mImageNew - mImage))))
%             disp(sprintf('alpha = %f', alpha));
            if (isinf(alpha))
                error('alpha is infinity');
            end
             % compute new alpha and new image and gradient
            alpha = alpha / beta;
            mImageNew = posProjection(mImage + alpha * mGrad);
        end
    end % end of if
    

end % end of function


function mMatNew = posProjection(mMat)
%
% Projects each entry of mMat to non-negative values.
%
    mMatNew = mMat;
    mMatNew(mMat < 0) = 0;
end



