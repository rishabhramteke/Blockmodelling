classdef HardMembershipIncAdj < HardMembershipInc
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

            % we do not need to acutally compute the frobenius norm, the L1 norm
            % will still tell us which position is minimal for a vertex to be
            % assigned to.
%             currObjVal = sum(sum(abs(mExistApprox)));


            % initialise the new membership matrix
            %mNewMembership = sparse(zeros(size(mMembership,1), size(mMembership,2)));

            % randomise order we evaluate the vertices
            vVertOrder = randperm(size(mMembership,1));

            % loop through each vertex/row and assign to best position
            for i = 1 : size(vVertOrder, 2)
            %     disp(sprintf('vertex %d', vVertOrder(i)));
        
                % store objective values
                vObjVals = zeros(1,size(mMembership,2));

                currV = vVertOrder(i);
                
                % reset positions for mMembership
                % mMembership(v,:) = zeros(1,size(mMembership,2));
                % find which position v belongs to currently
                currP = find(mMembership(currV,:) > 0);
                %     assert(~isempty(currP));
                % if singleton, we don't try moving
                if sum(mMembership(:, currP)) <= 1
                    continue;
                end

                
                for p = 1 : size(mMembership,2)
                    if (p == currP)
%                         vObjVals(p) = currObjVal;
                        vObjVals(p) = abs(sum(mExistApprox(:,currV)) + sum(mExistApprox(currV,:)) + mExistApprox(currV,currV));;
                    else
                        %mMembership(v,p) = 1;
                        %vObjVals(p) = objective(mAdj, mImage, mMembership);
                        %mMembership(v,p) = 0;
% %                         mNewApprox = mExistApprox + (mMI * mDeltaTrans + mDelta * mIM - mDelta * mImage * mDeltaTrans) .* obj.m_mWeight;
%                         mNewApprox = mExistApprox;
%                         mNewApprox(:,currV) = mNewApprox(:,currV) + (mMI(:,currP) - mMI(:, p)) .* obj.m_mWeight(:,currV);
%                         mNewApprox(currV,:) = mNewApprox(currV,:) + (mIM(currP,:) - mIM(p, :)) .* obj.m_mWeight(currV, :);
%                         mNewApprox(currV, currV) = mNewApprox(currV, currV) -...
%                             (mImage(currP, currP) - mImage(p, currP) - mImage(currP, p) + mImage(p,p)) * obj.m_mWeight(currV, currV);                        
%                         vObjVals(p) = sum(sum(abs(mNewApprox)));

                        diff = 0;
                        diff = diff + sum(mExistApprox(:,currV) + (mMI(:,currP) - mMI(:, p)) .* obj.m_mWeight(:,currV));
%                         diff = diff + sum(mExistApprox(:,currV) + obj.m_mWeight(:,currV)' * (mMI(:,currP) - mMI(:, p)));
%                         diff = diff + sum((mIM(currP,:) - mIM(p, :)) .* obj.m_mWeight(currV, :));
                        diff = diff + sum(mExistApprox(currV,:)' + (mIMTrans(:,currP) - mIMTrans(:,p)) .* mWeightTrans(:,currV));
                        diff = diff + mExistApprox(currV, currV) - (mImage(currP, currP) - mImage(p, currP) - mImage(currP, p) + mImage(p,p)) * obj.m_mWeight(currV, currV); 
                        vObjVals(p) = abs(diff);
                    end
                end % end of inner for
    
                % find minimum value
                minObjVal = min(vObjVals);
                vMinI = find(vObjVals == minObjVal);
                %     assert(size(vMinI,2) > 0);

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
                
%                 mMI = mMI - mDelta * mImage;
%                 mIM = mIM - mImage * mDeltaTrans;
                    mMI(currV,:) = mMI(currV,:) - (mImage(currP,:) - mImage(minP,:));
%                     mIM(:,currV) = mIM(:,currV) - (mImage(:,currP) - mImage(:,minP)); 
                    mIMTrans(currV,:) = mIMTrans(currV,:) - (mImageTrans(currP,:) - mImageTrans(minP,:));
                    % update objective value
%                     currObjVal = minObjVal;
                end
   
            end % end of outer for

        end % end of function
    
        
    end % end of methods
    
    
end % end of class



