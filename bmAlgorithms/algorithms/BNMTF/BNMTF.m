function [model, currObjVal, vObjVal, vImageDis, vKKTResidual, cvComparison, cmImage, cmMembership] =...
    BNMTF(G,loss_option,sparse_option,lambda, fImageDistanceFunc, cfValidMeasFunc, bDiscretiseMembership, max_rank, max_iter, varargin)
%G: the adjancency matrix of the network
%loss_option =0: square loss; 
%            =1: KL-divergence
%sparse_loss =0: all elements used; 
%            =1: positive elements used
%lambda: the regularization parameter
%max_rank: the maximum value for the rank (optional)
%U0: the initial guess for U (optional)
%B0: the initial guess for B (optional)
%model: struct variable with learned U and B as components
%
% Modified by: Jeffrey Chan, 2014
    


    % parse arguments
    inParser = inputParser;   
    inParser.KeepUnmatched = true;




    parse(inParser, varargin{:});
     


   
    % parse arguments
    inParser = inputParser;
    inParser.KeepUnmatched = true;

    % add parameters and set default valuse
    % initial membership and image matrices
    addParameter(inParser, 'initialMem', NaN);
    addParameter(inParser, 'initialImage', NaN);    
    % default is to collect per iteration stats
    addParameter(inParser, 'collectPerIterStats', false); 
    
    parse(inParser, varargin{:});
    
    U0 = inParser.Results.initialMem;
    B0 = inParser.Results.initialImage;
    bCollectPerIterStats = inParser.Results.collectPerIterStats;  
    

       

    N = size(G,1);
    
    if (nargin < 7)
        max_rank = ceil(N/2);
    end

    %remove non-connected nodes
    % Modified, do not remove non-connected nodes
    active_nodes = [1:size(G,1)];
%     active_nodes = sum(G)~=0;
    active_node_indices = find(active_nodes);
    non_active_nodes = ~active_nodes;
    non_active_node_indices = find(non_active_nodes);

    if sum(non_active_nodes)>0
        G = G(active_node_indices,active_node_indices);
        fprintf('\nNOTE: the following non-connected nodes are removed from the graphs:\n');
        disp(non_active_node_indices);
        if nargin>=6
            U0=U0(active_node_indices,:);
        end
    end
    
    if (nargin < 8)
        max_iter = 100;
    end
    
    %set up the initial W,H matrices
    if isnan(U0)
        U0 = rand(N,max_rank);
    end
    if isnan(B0)
        B0 = rand(max_rank);
        B0 = (B0+B0')/2;
    end    

    

    %run the BNMTF
%     max_iter=100;
    
    if loss_option==0
        fDistanceFunc = EucTriFact;
        if bCollectPerIterStats
            fKKTResidualFunc = EucKKTResidual(fDistanceFunc);      
            [U1, B, vObjVal, vImageDis, vKKTResidual, cvComparison, cmImage, cmMembership] =...
                BNMTF_sq(G, U0, B0, lambda, sparse_option, max_iter, fDistanceFunc, fImageDistanceFunc, cfValidMeasFunc, 'kktFunc', fKKTResidualFunc, 'innerCollectPerIterStats', true);            
        else
            [U1, B, vObjVal, vImageDis, vKKTResidual, cvComparison, cmImage, cmMembership] =...
                BNMTF_sq(G, U0, B0, lambda, sparse_option, max_iter, fDistanceFunc, fImageDistanceFunc, cfValidMeasFunc, 'innerCollectPerIterStats', false);
        end

        currObjVal = fDistanceFunc.distance(G, B, U1);
    elseif loss_option==1
        [U1, B] = BNMTF_kl(G, U0, B0, lambda, sparse_option, max_iter, 0.1);
    else
        error('Wrong loss option!');
    end
    
    if sum(non_active_nodes)>0
        U=zeros(N,max_rank);
        U(active_node_indices,:)=U1;
    else
        U=U1;
    end
    
    % see if we need to discretise the results
    if bDiscretiseMembership
       U = discretise(U);
       % loss options should have passed and these objects should have been
       % created
       % compute the new distances
       currObjVal = fDistanceFunc.distance(G, B, U);
       vObjVal(end) = currObjVal;
       if bCollectPerIterStats
            vKKTResidual(end) = fKKTResidualFunc.residual(G, B, U);
       end
       cmMembership{end} = U;       
    end        
    
    model.U=U;
    model.B=B;
    clear U B U1 B1 U0 B0 active_nodes active_node_indices non_active_nodes non_active_node_indices;    