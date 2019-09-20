classdef SoftMembershipCoordAdj < SoftMembershipCoord
    %
    % Soft membership coordinate update.
    % Based on "Overlapping Community Detection via Bounded nonnegative matrix
    % Tri-factorisation"
    %
    %
    % @author: Jeffrey Chan, 2013
    %
    
    properties
        % pre-computed weight matrix.
        m_mWeight;
    end
   
    methods
    
        function obj = SoftMembershipCoordAdj(mWeight)
            obj.m_mWeight = mWeight;
        end
       
        
        function [mMembership] = updateMembership(obj, mAdj, mImage, mMembership)
            %
            % Updates the membership according to the coordinate descent
            % approach described in 'overlapping community detection via bounded
            % nonnegative matrix factorization'
            %
            

            % Code starting from here is the work of BNMTF
            %%%%%%%
            [n,k]=size(mMembership);
            U=mMembership;
            B=mImage;
            
            tmp1=full((U*B*U'-mAdj) .* obj.m_mWeight) ;
            tmp2=U*B;
            tmp3=B*U';

            %%%Update U
            for p=1:n
                for q=1:k
                    %%%Find t
                    Z1=sparse(n,n);
                    Z1(:,p)=tmp2(:,q);
                    Z1(p,:)= Z1(p,:)+tmp3(q,:);
                    Z1 = Z1 .* obj.m_mWeight;
                    
                    Z1_index=find(Z1>0);
                    z2=B(q,q) * obj.m_mWeight(p,q);
                    
                    a=z2^2;
                    b=full(2*z2*Z1(p,p)); 
                    c=full(2*z2*tmp1(p,p)+norm(Z1,'fro')^2);
                    d=full(2*tmp1(Z1_index)'*Z1(Z1_index));
                    
                    [t,flag]= obj.find_optimal_U(a,b,c,d,[-U(p,q),1-U(p,q)]);
                    clear Z1_index a b c d;
                    %%%Update tmp1, tmp2, tmp3
                    if t~=0
                        tmp1=full(tmp1+t*Z1);
                        tmp1(p,p)=tmp1(p,p)+(t^2)*z2;
                        tmp2(p,:)=tmp2(p,:)+t*B(q,:);
                        tmp3(:,p)=tmp3(:,p)+t*B(:,q);
                    end
                    clear Z1 z2;
                    %%%Update U(p,q)
                    if flag==-1
                        U(p,q)=U(p,q)+t;
                    elseif flag==0
                        if abs(U(p,q)+t)>sqrt(eps)
                            error('wrong left flag');
                        end
                        U(p,q)=0;
                    else
                        if abs(U(p,q)+t-1)>sqrt(eps)
                            error('wrong right flag');
                        end
                        U(p,q)=1;
                    end
                end
            end
            
            %%%%%%%%%%
    
        end % end of function
        
        
        
    
        
    end % end of methods
    
    
end % end of class




