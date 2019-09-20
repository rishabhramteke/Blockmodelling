classdef EucTriFactBM < EucTriFact
    % 
    % Euclidean tri-factorisation with sigmiod constraint.
    %
    %
    % @author: Jeffrey Chan, 2013
    %    
    
    % sigmiod parameters
    properties (GetAccess = 'public', SetAccess = 'protected')
        % default settings
%         m_gamma = 0.25;
%         m_beta_mu = 400;
%         m_mTau = 0.5;
        % function class to compute the various regularisation terms
        % Note this default uses the default parameters for SigmoidImageReg.
        % Initialise with a function object if customised parameters.
        m_fImageReg = SigmoidImageReg(0.001, 1000);
    end % end of properties
    
    % regularisation parameter
    properties 
        m_regPara = 100;
    end % end of properties
    
    % member objects
    properties
%         m_objEuc = EucTriFact;
    end

    
    methods
        function obj = EucTriFactBM(varargin)
            %
            % constructor
            %
            
            
            obj@EucTriFact(varargin{:});
            
            inParser = inputParser;
            inParser.KeepUnmatched = true;
            
            addParameter(inParser, 'imageRegFunc', SigmoidImageReg);
            addParameter(inParser, 'imageRegWeight', 100);

            parse(inParser, varargin{:});            
            
            obj.m_fImageReg = inParser.Results.imageRegFunc;            
            obj.m_regPara = inParser.Results.imageRegWeight;   

%             if ~isempty(varargin)
%                 if length(varargin) >= 1
%                     obj.m_fImageReg = varargin{1}; 
%                 end                
%                 if length(varargin) >= 2
%                     obj.m_regPara = varargin{2}; 
%                 end
% %                 if length(varargin) >= 3
% %                     obj.m_beta_mu = varargin{3};
% %                 end
% %                 if length(varargin) == 4
% %                     obj.m_mTau = varargin{4};
% %                 end
%             end    
            
            
        end
        
        
%         function adjustIdealShift(obj, mTau)
%             obj.m_mTau = mTau;
%         end
        
        
        function [dis] = distance(obj, mAdj, mImage, mMembership, varargin)
            %
            % Compute the approximation distance between mAdj and the
            % approximated blockmodel + regularisation term.
            %

            % rescaling factor for regularisation term
%             rescalFactor = (size(mAdj,1)^2 / size(mImage,1)^2);
            
%             mIdeal = 1 ./ (1 + obj.m_gamma * exp(-obj.m_beta_mu * (mImage - obj.m_tau)));
%             mIdeal = obj.computeIdeal(mImage);
%             mApprox = mIdeal - mImage;
%             objVal = trace(mApprox * mApprox');
            objVal =  obj.m_fImageReg.computeImageRegObj(mImage);
    
            dis = obj.m_regPara *  objVal;
            
            % see if there is a weight matrix
            if length(varargin) == 1
                mWeight = varargin{1};
                dis = dis + distance@EucTriFact(obj, mAdj, mImage, mMembership, mWeight);
%                 dis = dis + obj.m_objEuc.distance(mAdj, mImage, mMembership, mWeight);
            else
                dis = dis + distance@EucTriFact(obj, mAdj, mImage, mMembership);
%                 dis = dis + obj.m_objEuc.distance(mAdj, mImage, mMembership);
            end
            
            % convert to number if in sparse representation.
            if issparse(dis)
                dis = full(dis);
            end
        end % end of distance()
        
        
        function dis = imageRegObjective(obj, mImage)
            dis = obj.m_fImageRef.computeImageRegObj(mImage);
        end % end of function
        
%         
%         function mIdeal = computeIdeal(obj, mImage)
%             %
%             % Computes the ideal matrix (sigmoid) function.
%             %
%             
%             if obj.m_sImageReg == 'gaussian'
%                 mIdeal = exp(-0.5 * (mImage - obj.m_gaussMean)
%             else
%                 mIdeal = 1 ./ (1 + obj.m_gamma * exp(-obj.m_beta_mu * (mImage - obj.m_mTau)));
%         end
        
        
        
        function [mGrad] = gradient(obj, mAdj, mImage, mMembership)
            %
            % Compute the gradient of approximation function wrt to mImage.
            %
            
            % rescaling factor for regularisation term
%             rescalFactor = (size(mAdj,1)^2 / size(mImage,1)^2);
            
            mFirst = -2 * mMembership' * mAdj * mMembership;
            mMemTerm = mMembership' * mMembership;
            mSecond = 2 * (mMemTerm * mImage * mMemTerm);
            
%             mRegGrad = obj.regularisationGradient(mImage);
            mRegGrad = obj.imageRegGradient(mImage);
            
            % add all the pieces intogether
            mGrad = mFirst + mSecond + obj.m_regPara * mRegGrad;
        end % end of function gradient()
        
        
        function [mRegGrad] = imageRegGradient(obj, mImage)
            %
            % Compute the gradient of the image regularisation term.
            % Does not factor in the regularisation weight.
            %
%             
%             mExp = obj.m_gamma * exp(-obj.m_beta_mu *(mImage - obj.m_mTau));
%             mExpSq = (1 + mExp) .* (1 + mExp);
%             mThird = 1 + mExp .* (1 + obj.m_beta_mu * mImage) ./ mExpSq;
%             
%             
%             mExpCb = mExpSq .* (1 + mExp);
%             mFifth = 2 * obj.m_beta_mu * (mExp ./ mExpCb);
%             
%             mRegGrad = mThird + 2 * mImage + mFifth;

            mRegGrad = obj.m_fImageReg.computeImageRegGrad(mImage);
        end % end of function()
        

        function [dis] = computeImageRegDis(obj, mAdj, mImage, stepSize, basisRow, basisCol)
%             rescalFactor = (size(mAdj,1)^2 / size(mImage,1)^2);
            dis = obj.m_regPara  * obj.imageRegCoordDesc(mImage, stepSize, basisRow, basisCol);
        end
        
        
        function [dis] = coordDescDistance(obj, mAdj, mImage, mMembership, stepSize, basisRow, basisCol)
            %
            % Compute the objective value of the coordinate descent expression.
            % TODO: we can reuse most of the calculations again, need to figure
            % out how to store them between calls.
            %
            
            % rescaling factor for regularisation term
%             rescalFactor = (size(mAdj,1)^2 / size(mImage,1)^2);
            
            
            % construct the unit basis
%             mUnitBasis = zeros(size(mImage,1), size(mImage,2));
%             mUnitBasis(basisRow, basisCol) = 1;
%             mUnitBasis = sparse(mUnitBasis);
%             
%             mX0 = mAdj - mMembership * mImage * mMembership';
%             mX1 = mMembership * mUnitBasis * mMembership';
%             disSq = stepSize * stepSize;
%             % Tr(X0' X0)
%             dis = trace(mX0' * mX0);
%             % Tr(X0' X1)
%             dis = dis - 2 * stepSize * trace(mX0' * mX1);
%             % Tr(X1' X1)
%             dis = dis + disSq * trace(mX1' * mX1);
            dis = coordDescDistance@EucTriFact(obj, mAdj, mImage, mMembership, stepSize, basisRow, basisCol);
            
%             mX2 = mImage - obj.m_mTau;
%             % compute the term where we don't add t
%             mExp = 1 ./ (1 + obj.m_gamma * exp(-obj.m_beta_mu * mX2));
%             % subtract the exponential term that isn't included in the sum
%             mExp(basisCol, basisRow) = 0;
%             
%             % compute the exp term where we add stepSize
%             singleExp = 1 ./ (1 + obj.m_gamma * exp(-obj.m_beta_mu * (stepSize + mX2(basisRow, basisCol))));
            
            % compute the whole sum
%             dis = dis + obj.m_regPara * rescalFactor * (disSq + singleExp * singleExp - 2 * stepSize * singleExp + 2 * stepSize * mImage(basisRow, basisCol)...
%                 - 2 * mImage(basisRow, basisCol) * singleExp + trace(mExp' * (mExp + mImage)) + trace(mImage' * mImage));
            dis = dis + obj.m_regPara * obj.imageRegCoordDesc(mImage, stepSize, basisRow, basisCol);
        end % end of function coordDescDistance
        
        
        function dis = imageRegCoordDesc(obj, mImage, stepSize, basisRow, basisCol)
            dis = obj.m_fImageReg.computeImageRegCoord(mImage, stepSize, basisRow, basisCol);
        end % end of function
        
    
    end




end % end of class