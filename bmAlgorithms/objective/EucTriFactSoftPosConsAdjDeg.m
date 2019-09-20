classdef EucTriFactSoftPosConsAdjDeg < EucTriFactSofPosConsAdj
    % 
    % Euclidean tri-factorisation with adjusted
    % edge mismatch penalties.
    %
    % Degree null model.
    %
    %
    % @author: Jeffrey Chan, 2013
    %    
    
    
    
    properties 
        % in degree of vertices
        m_vInDeg;
        % out degree of vertices
        m_vOutDeg;
        % in and out degree for each block
        m_vBlockInDeg;
        m_vBlockOutDeg;
    end % end of properties

    
    
    methods
        function obj = EucTriFactSoftPosConsAdjDeg(mWeight, varargin)
            %
            % constructor
            % mWeight - weight adjustment matrix.
            %
            
            obj@EucTriFactSofPosConsAdj(mWeight, varargin{:});
           
        end
        
        
        function initDegrees(obj, mAdj, mMembership)
            %
            % Initialise the degrees.
            %
            
            for v = 1 : size(mAdj,1)
               obj.m_vInDeg(v) = sum(mAdj(v,:));
               obj.m_vOutDeg(v) = sum(mAdj(:,v));
            end
            
            obj.m_vBlockInDeg = mMembership' * obj.m_vInDeg;
            obj.m_vBlockOutDeg = mMembership' * obj.m_vOutDeg;
        end
        
        
        function updateDegree(obj, mMembership, v, origPos, membershipDecr, newPos, membershipIncr)
            %
            % Move vertex v from origPos (decrease in membership) to newPos
            % (increase in membership).
            %
            % Function assumes that checks have been performed for the
            % membershipDecr and membershipIncr.
            %
            
            
            obj.m_vBlockInDeg(origPos) = obj.m_vBlockInDeg(origPos) - membershipDecr * obj.m_vInDeg(v);
            obj.m_vBlockOutDeg(origPos) = obj.m_vBlockOutDeg(origPos) - membershipDecr * obj.m_vOutDeg(v);
            
            obj.m_vBlockInDeg(newPos) = obj.m_vBlockInDeg(newPos) + membershipIncr * obj.m_vInDeg(v);
            obj.m_vBlockOutDeg(newPos) = obj.m_vBlockOutDeg(newPos) + membershipIncr * obj.m_vOutDeg(v);
            
            % update out degree
            for p = 1 : size(mMembership, 2)    
                obj.m_mWeight(origPos, p) = obj.m_vBlockOutDeg(origPos) * obj.m_vBlockInDeg(p);
                obj.m_mWeight(newPos, p) = obj.m_vBlockOutDeg(newPos) * obj.m_vBlockInDeg(p);
            end
            
            % update in degree
            for p = 1 : size(mMembership, 2)    
                if p == origPos || p == newPos
                    continue;
                end
                obj.m_mWeight(p, origPos) = obj.m_vBlockOutDeg(p) * obj.m_vBlockInDeg(origPos);
                obj.m_mWeight(p, newPos) = obj.m_vBlockOutDeg(p) * obj.m_vBlockInDeg(newPos);
            end
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function [dis] = distance(obj, mAdj, mImage, mMembership)
            %
            % Compute the approximation distance between mAdj and the
            % approximated blockmodel + regularisation term.
            %

            dis = distance@EucTriFactSofPosConsAdj(obj, mAdj, mImage, mMembership, obj.m_mWeight);
        end % end of distance()
        
        
        
        function [mGrad] = gradient(obj, mAdj, mImage, mMembership)
            %
            % Compute the gradient of approximation function.
            %
            

            
            % this is the term that is different from EucTriFactBM
            mGrad = 2 * mMembership' * ((mMembership * mImage * mMembership' - mAdj) .* obj.m_mWeight.^2) * mMembership;
        end % end of function gradient()
        
        
        
        
        function [dis] = coordDescDistance(obj, mAdj, mImage, mMembership, stepSize, basisRow, basisCol)
            %
            % Compute the objective value of the coordinate descent expression.
            % TODO: we can reuse most of the calculations again, need to figure
            % out how to store them between calls.
            %

            
            % construct the unit basis
            mUnitBasis = zeros(size(mImage,1), size(mImage,2));
            mUnitBasis(basisRow, basisCol) = 1;
            mUnitBasis = sparse(mUnitBasis);
            
            mX0 = (mAdj - mMembership * mImage * mMembership') .* obj.m_mWeight;
            mX1 = (mMembership * mUnitBasis * mMembership') .* obj.m_mWeight;
            disSq = stepSize * stepSize;
            % Tr(X0' X0)
            dis = trace(mX0' * mX0);
            % Tr(X0' X1)
            dis = dis - 2 * stepSize * trace(mX0' * mX1);
            % Tr(X1' X1)
            dis = dis + disSq * trace(mX1' * mX1);
        end % end of function coordDescDistance
        
    
    end




end % end of class