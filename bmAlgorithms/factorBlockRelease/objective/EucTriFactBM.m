classdef EucTriFactBM < EucTriFact
    % 
    % Euclidean tri-factorisation with regularisation term.
    %
    % If using this code, please kindly cite the following paper:
    % J. Chan, W. Liu, C. Leckie, J. Bailey and K. Ramamohanarao. 
    % "Discovering Latent Blockmodels in Sparse and Noisy Graphs using Non-Negative Matrix Factorisation."
    % In Proceedings of 22nd ACM International Conference on Information and Knowledge Management, October 2013.            
    %
    %
    % @author: Jeffrey Chan, 2013
    %    
    
    % sigmiod parameters
    properties (GetAccess = 'public', SetAccess = 'protected')
        % function class to compute the various regularisation terms
        % Note this default uses the default parameters for SigmoidImageReg.
        % Initialise with a function object if customised parameters.
        m_fImageReg = SigmoidImageReg;
    end % end of properties
    
    % regularisation parameter
    properties 
        m_regPara = 100;
    end % end of properties
    
    % member objects
    properties
        m_objEuc = EucTriFact;
    end

    
    methods
        function obj = EucTriFactBM(varargin)
            %
            % constructor
            %
            
            if ~isempty(varargin)
                if length(varargin) >= 1
                    obj.m_fImageReg = varargin{1}; 
                end                
                if length(varargin) >= 2
                    obj.m_regPara = varargin{2}; 
                end
            end    

        end
        
        
%         function adjustIdealShift(obj, mTau)
%             obj.m_mTau = mTau;
%         end
        
        
        function [dis] = distance(obj, mAdj, mImage, mMembership, varargin)
            %
            % Compute the approximation distance between mAdj and the
            % approximated blockmodel + regularisation term.
            %

            objVal =  obj.m_fImageReg.computeImageRegObj(mImage);
    
            dis = obj.m_regPara *  objVal;
            
            % see if there is a weight matrix
            if length(varargin) == 1
                mWeight = varargin{1};
                dis = dis + obj.m_objEuc.distance(mAdj, mImage, mMembership, mWeight);
            else
                dis = dis + obj.m_objEuc.distance(mAdj, mImage, mMembership);
            end
            
            % convert to number if in sparse representation.
            if issparse(dis)
                dis = full(dis);
            end
        end % end of distance()
        
        
        function dis = imageRegObjective(obj, mImage)
            dis = obj.m_fImageRef.computeImageRegObj(mImage);
        end % end of function
        

               
        
        function [mGrad] = gradient(obj, mAdj, mImage, mMembership)
            %
            % Compute the gradient of approximation function.
            %
            

            
            mFirst = -2 * mMembership' * mAdj * mMembership;
            mMemTerm = mMembership' * mMembership;
            mSecond = 2 * (mMemTerm * mImage * mMemTerm);
            
            mRegGrad = obj.imageRegGradient(mImage);
            
            % add all the pieces intogether
            mGrad = mFirst + mSecond + obj.m_regPara * mRegGrad;
        end % end of function gradient()
        
        
        function [mRegGrad] = imageRegGradient(obj, mImage)
            %
            % Compute the gradient of the image regularisation term.
            % Does not factor in the regularisation weight.
            %           


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
            

            dis = coordDescDistance@EucTriFact(obj, mAdj, mImage, mMembership, stepSize, basisRow, basisCol);
            
            dis = dis + obj.m_regPara * obj.imageRegCoordDesc(mImage, stepSize, basisRow, basisCol);
        end % end of function coordDescDistance
        
        
        function dis = imageRegCoordDesc(obj, mImage, stepSize, basisRow, basisCol)
            dis = obj.m_fImageReg.computeImageRegCoord(mImage, stepSize, basisRow, basisCol);
        end % end of function
        
    
    end




end % end of class