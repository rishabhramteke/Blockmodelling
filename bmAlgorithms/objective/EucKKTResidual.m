classdef EucKKTResidual
    %
    % Class of function object used to calculate the KKT residual of the
    % Euclidean loss functions.
    %
    % @author: Jeffrey Chan, 2014
    %    
    
    properties
       m_fEucTriFact;
    end
    

    methods
        function obj = EucKKTResidual(fEucDisFunc)
            obj.m_fEucTriFact = fEucDisFunc;
            
        end
        
        
        function [residualVal] = residual(obj, mAdj, mImage, mMembership)
            %
            % Compute the KKT residual using the Euclidean norm.
            %
            
            % compute the gradient of mMembership and mImage
            mPosGrad = obj.m_fEucTriFact.posGradient(mAdj, mImage, mMembership);
            mImageGrad = obj.m_fEucTriFact.gradient(mAdj, mImage, mMembership);
            
            
            if or(~isreal(mMembership), isnan(mMembership))
                display(mMembership);
            end
            
            if or(~isreal(mImage), isnan(mImage))
                display(mImage);
            end
            
            if or(~isreal(mPosGrad), isnan(mPosGrad))
                display(mPosGrad);
            end
            
            if or(~isreal(mImageGrad), isnan(mImageGrad))
                display(mImageGrad);
            end
            
            mProjPosGrad = ((mPosGrad < 0) | mMembership) .* mPosGrad;
            mProjImageGrad = ((mImageGrad < 0) | mImage) .* mImageGrad;
            
            residualVal = sqrt(norm(mProjPosGrad,2)^2 + norm(mProjImageGrad, 2)^2);
        end % end of function
        
    end % end of methods

end % end of class
