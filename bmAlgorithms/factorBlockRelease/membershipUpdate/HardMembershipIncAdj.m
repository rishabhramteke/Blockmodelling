classdef HardMembershipIncAdj < HardMembershipInc
    %
    % Incremental hard membership update, adjusted for sparsity.
    %
    % If using this code, please kindly cite the following paper:
    % J. Chan, W. Liu, C. Leckie, J. Bailey and K. Ramamohanarao. 
    % "Discovering Latent Blockmodels in Sparse and Noisy Graphs using Non-Negative Matrix Factorisation."
    % In Proceedings of 22nd ACM International Conference on Information and Knowledge Management, October 2013.
    %
    % @author: Jeffrey Chan, 2013
    %
    
    properties
        % pre-computed weight matrix.
        m_mWeight;
    end
   
    methods
    
        function obj = HardMembershipIncAdj(mWeight)
            obj.m_mWeight = mWeight;
        end
        
        
        function [mMembership] = updateMembership(obj, mAdj, mImage, mMembership)
            %
            % Updates the hard membership, using a greedy, vertex by vertex assignment.
            % Ordering of update does not affect final result.
            %
            % Incremental update.
            %
            % This version will multiple by an adjustment weight.


            % precompute some matrices
            mMI = mMembership * mImage;
            % we work on the transpose of mImage * mMembership', for speed
            mImageTrans = mImage';
            mIMTrans = mMembership * mImageTrans;
            mExistApprox = (mAdj - mMI * mMembership') .* obj.m_mWeight;
            mWeightTrans = obj.m_mWeight';



            % randomise order we evaluate the vertices
            vVertOrder = randperm(size(mMembership,1));

            % loop through each vertex/row and assign to best position
            for i = 1 : size(vVertOrder, 2)
        
                % store objective values
                vObjVals = zeros(1,size(mMembership,2));

                currV = vVertOrder(i);
                
                % find which position v belongs to currently
                currP = find(mMembership(currV,:) > 0);
                % if singleton, we don't try moving
                if sum(mMembership(:, currP)) <= 1
                    continue;
                end

                
                for p = 1 : size(mMembership,2)
                    if (p == currP)
                        vObjVals(p) = abs(sum(mExistApprox(:,currV)) + sum(mExistApprox(currV,:)) + mExistApprox(currV,currV));
                    else

                        diff = 0;
                        diff = diff + sum(mExistApprox(:,currV) + (mMI(:,currP) - mMI(:, p)) .* obj.m_mWeight(:,currV));
                        diff = diff + sum(mExistApprox(currV,:)' + (mIMTrans(:,currP) - mIMTrans(:,p)) .* mWeightTrans(:,currV));
                        diff = diff + mExistApprox(currV, currV) - (mImage(currP, currP) - mImage(p, currP) - mImage(currP, p) + mImage(p,p)) * obj.m_mWeight(currV, currV); 
                        vObjVals(p) = abs(trace(diff));
                    end
                end % end of inner for
    
                % find minimum value
                minObjVal = min(vObjVals);
                vMinI = find(vObjVals == minObjVal);
                
                % only change assignment if minI != currP 
                if isempty(find(vMinI == currP, 1))
                    % pick one of the minimum
                    % randsample can be slow to continually call, so we only
                    % call it when we have to
                    if length(vMinI) > 2
                        minP = vMinI(randsample(length(vMinI), 1));
                    else   
                        minP = vMinI(1);
                    end                    
                   
                    % update mMembership
                    mMembership(currV, currP) = 0;
                    mMembership(currV, minP) = 1;
                    % update all the precalculated matrices                    
                    mExistApprox(:,currV) = mExistApprox(:,currV) + (mMI(:,currP) - mMI(:, minP)) .* obj.m_mWeight(:,currV);
                    mExistApprox(currV,:) = mExistApprox(currV,:) + ((mIMTrans(:,currP) - mIMTrans(:,p)) .* mWeightTrans(:,currV))';
                    mExistApprox(currV, currV) = mExistApprox(currV, currV) -...
                        (mImage(currP, currP) - mImage(minP, currP) - mImage(currP, minP) + mImage(minP,minP)) * obj.m_mWeight(currV, currV);                
                
                    mMI(currV,:) = mMI(currV,:) - (mImage(currP,:) - mImage(minP,:));
                    mIMTrans(currV,:) = mIMTrans(currV,:) - (mImageTrans(currP,:) - mImageTrans(minP,:));
                end
   
            end % end of outer for

        end % end of function
    
        
    end % end of methods
    
    
end % end of class



