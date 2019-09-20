classdef HardMembershipInc
    %
    % @author: Jeffrey Chan, 2013
    %
   
methods
    
    function [mMembership] = updateMembership(obj, mAdj, mImage, mMembership)
    %
    % Updates the hard membership, using a greedy, vertex by vertex assignment.
    % Ordering of update does not affect final result.
    %
    % Incremental update.
    %



        % precompute some matrices
        mMI = mMembership * mImage;
%         mIM = mImage * mMembership';
        mImageTrans = mImage';
        mIMTrans = mMembership * mImageTrans;        
        mExistApprox = mAdj - mMI * mMembership';

        % we do not need to compute the multiplication and trace, as
        % we only want the minimal one, so we can keep in L1 norm
%         currObjVal = trace(mExistApprox * mExistApprox');
%         currObjVal = sum(sum(abs(mExistApprox)));

        % initialise two sparce vectors to quickly update memberships
%         mDelta = sparse(zeros(size(mMembership,1), size(mMembership,2)));

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
            %disp(currP);
            %     assert(~isempty(currP));
            % if singleton, we don't try moving
            if sum(mMembership(:, currP)) <= 1
                continue;
            end

            
%             mDelta(vVertOrder(i), currP) = 1;
            
    
            for p = 1 : size(mMembership,2)
                disp(currP);
                disp(p);
                if (p == currP)
                %if (p == 2)
                    disp('if');
%                     vObjVals(p) = currObjVal;
                    vObjVals(p) = abs(sum(mExistApprox(:,currV)) + sum(mExistApprox(currV,:)) + mExistApprox(currV,currV));
                else
                     disp('else');
                    %mMembership(v,p) = 1;
                    %vObjVals(p) = objective(mAdj, mImage, mMembership);
                    %mMembership(v,p) = 0;
%                     mDelta(currV, p) = -1;
%                     mDeltaTrans = mDelta';
%                     mNewApprox = mExistApprox + mMI * mDeltaTrans + mDelta * mIM - mDelta * mImage * mDeltaTrans;
%                     mNewApprox = mExistApprox + mMI * mDeltaTrans + mDelta * mIM;
%                     mNewApprox = mExistApprox + mMI * mDeltaTrans;
%                     mNewApprox = mExistApprox;
%                     mNewApprox(:,currV) = mNewApprox(:,currV) + (mMI(:,currP) - mMI(:, p));
%                     mNewApprox(currV,:) = mNewApprox(currV,:) + (mIM(currP,:) - mIM(p, :));
%                     mNewApprox(currV, currV) = mNewApprox(currV, currV) -...
%                         (mImage(currP, currP) - mImage(p, currP) - mImage(currP, p) + mImage(p,p));
                    % we do not need to compute the multiplication and trace, as
                    % we only want the minimal one, so we can keep in L1 norm
                    % (but have to take absolute value).
%                     vObjVals(p) = trace(mNewApprox * mNewApprox');
%                     vObjVals(p) = sum(sum(abs(mNewApprox)));
                    
                   diff = sum(mExistApprox(:,currV) + mMI(:,currP) - mMI(:, p));
                   diff = diff + sum(mExistApprox(currV,:)' + mIMTrans(:,currP) - mIMTrans(:,p));
                   diff = diff + mExistApprox(currV, currV) - (mImage(currP, currP) - mImage(p, currP) - mImage(currP, p) + mImage(p,p)); 
                   % we do not need to compute the L2 norm of the part that changes, as
                   % we only want the minimal one, so we can keep in L1 norm
                   % (but have to take absolute value).
                   disp(abs(diff));
                   disp(vObjVals(p));
                   disp(p);
                   vObjVals(p) = abs(trace(diff));                    
                    
%                     mDelta(currV, p) = 0;
                end
            end % end of inner for
%             mDelta(currV, currP) = 0;
    
            % find minimum value
            minObjVal = min(vObjVals);
            vMinI = find(vObjVals == minObjVal);
%             if size(vMinI,2) <= 0
%                 vObjVals
%                 minObjVal
%                 vMinI
%             end

            % only change assignment if minI != currP 
            if isempty(find(vMinI == currP, 1))
                % pick one of the minimum
                % randsample can be slow to continually call, so we only
                % call it when we have to
                if length(vMinI) >= 2
                    minP = vMinI(randsample(length(vMinI), 1));
                else   
                    minP = vMinI(1);
                end
                
%                 minP = vMinI(randsample(size(vMinI,2), 1));
                % update mMembership
                mMembership(currV, currP) = 0;
                mMembership(currV, minP) = 1;
                % update all the precalculated matrices
%                 
%                 mDelta(currV, minP) = -1;
%                 mDeltaTrans = mDelta';
%                 mExistApprox = mExistApprox + mMI * mDeltaTrans + mDelta * mIM - mDelta * mImage * mDeltaTrans;
                mExistApprox(:,currV) = mExistApprox(:,currV) + (mMI(:,currP) - mMI(:, minP));
                mExistApprox(currV,:) = mExistApprox(currV,:) + (mIMTrans(:,currP) - mIMTrans(:,p))';
                mExistApprox(currV, currV) = mExistApprox(currV, currV) -...
                        (mImage(currP, currP) - mImage(minP, currP) - mImage(currP, minP) + mImage(minP,minP));                
                
%                 mMI = mMI - mDelta * mImage;
%                 mIM = mIM - mImage * mDeltaTrans;
%                 mMI(currV,:) = mMI(currV,:) - (mImage(currP,:) - mImage(minP,:));
%                 mIM(:,currV) = mIM(:,currV) - (mImage(:,currP) - mImage(:,minP));
                mMI(currV,:) = mMI(currV,:) - (mImage(currP,:) - mImage(minP,:));
                mIMTrans(currV,:) = mIMTrans(currV,:) - (mImageTrans(currP,:) - mImageTrans(minP,:));
                % update objective value
%                 currObjVal = minObjVal;
            end
   
        end % end of outer for

    end % end of function
    
    
        
end % end of methods
    
    
end % end of class



