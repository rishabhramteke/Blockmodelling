function [mImage] = updateCoordDescImageExact(mAdj, mImage, mMembership, fObjective)
%
% Updates mImage according to a coordinate descent update rule, with exact step
% size calculated (no need for line search).
% Should only be used for the euclidean distance.
%
%
% @author: Jeffrey Chan, 2013
%
    
%     mUnitBasis = sparse(zeros(size(mImage,1), size(mImage,2)));

    % euclidean
    mX0 = mAdj - mMembership * mImage * mMembership';
    vX0Trans = mX0(:)';
    
    % loop through coordinates
    for r = 1 : size(mImage,1)
        for c = 1 : size(mImage,2)
%             mUnitBasis(r,c) = 1;
            %mX1 = mMembership * mUnitBasis * mMembership';
            % outer product
            mX1 = mMembership(:,r) * mMembership(:,c)';
%             vX1 = mX1(:);
            % computing trace(mX1' * mX1)
            denom = full(sum(sum(mX1.^2)));
%             denom = full(vX1' * vX1);
%             fprintf('denom = %.2f\n', denom);
%             fprintf('trace() = %.2f\n', trace(mX1' * mX1));
            
            nomin = full(vX0Trans * mX1(:));
%             fprintf('nomin = %.2f\n', nomin);
%             fprintf('trace() = %.2f\n', trace(mX0' * mX1));            
            % compute optimal step size
%             altStepSize = trace(mX0' * mX1)/trace(mX1' * mX1);
            altStepSize = nomin / denom;
            minStepSize = min(1 - mImage(r,c), altStepSize);
            % update mImage
            mImage(r,c) = minStepSize + mImage(r,c);
%             mUnitBasis(r,c) = 0;
        end % end of inner for
        
    end % end of outer for
    

end % end of function











